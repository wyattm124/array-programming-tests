#include <array>
#include <complex>
#include <cstddef>
#include <cmath>
#include <arm_neon.h>
#include "prime_factor.hpp"
#include <iostream>

/// TODO: 
/// (1) - Raders DFT implementations for 3 - 7
/// (2) - specialized DFT implementations for 2, 4, 6, 8
/// (3) - implement these as base cases in FFTPlan
/// (4) - implemnt them for recursive steps in FFT plan
/// (5) - implement Bluestiens's for numbers that do not decompse into these factors
/// (6) - implement a top level search for best number to pad an input to for FFT

// For LLVM MCA Analysis
#define MCA_START __asm volatile("# LLVM-MCA-BEGIN");
#define MCA_END __asm volatile("# LLVM-MCA-END");

namespace FFT {
    // Neon intrinsics
    inline float32x2_t neon_mul(const float32x2_t &__restrict__ a, const float32x2_t &__restrict__ b) noexcept {
        return {a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
    }
    inline float neon_abs(const float32x2_t &__restrict__ a) noexcept {
        return std::sqrt(a[0]*a[0] + a[1]*a[1]);
    }
    inline float32x2_t neon_conj(const float32x2_t &__restrict__ a) noexcept {
        return {a[0], -a[1]};
    }

    // For DFT Transpositions
    template <unsigned int A, unsigned int B, typename T>
    void DFT_binner(T *__restrict__ in, T *__restrict__ out) noexcept {
        // Base case
        if constexpr (A * B == 1) {
            return;
        }

        // This will create B DFTs of size A
        for (std::size_t n = 0; n < A * B; n++) {
            const auto i = n % B;
            const auto j = n / B;
            const auto k = i * A + j;
            out[k] = in[n];
        }
        return;
    }
    
    // NOTE : although this algorithm is in-place, all the added computation is not worth
    //  the memory savings.

    // This will create B DFTs of size A
    template <unsigned int A, unsigned int B>
    void DFT_binner_inplace(float32x2_t *__restrict__ data) {
        // Base case
        if constexpr (A < 2) {
            return;
        }

        // Have an array to track which elements have already been placed
        std::array<bool, (A * B)> flip_tracker = {false};

        // Given an index n, return the index it should be moved to
        constexpr auto get_next_index = [=](std::size_t n) -> std::size_t {
            const auto i = n % A;
            const auto j = n / A;
            return i * B + j;
        };

        // Because of basic group theory, if N / p is not prime, we will
        //  have at least one loop of transposition that does not span all items.
        //  so as we whack-a-mole replace elements, we may need to find a new
        //  loop of elements to start replacing.
        std::size_t earliest_unmoved_index = 0;
        std::size_t curr_index = 0;
        std::size_t next_index = get_next_index(curr_index);
        float32x2_t temp_val = data[curr_index];
        while (earliest_unmoved_index < (A * B)) {
            // Complete one loop of transpositions
            while (!flip_tracker[next_index]) {
                flip_tracker[next_index] = true;
                if (curr_index != next_index) { // avoid swap with self
                    std::swap(temp_val, data[next_index]);
                    curr_index = next_index;
                    next_index = get_next_index(curr_index);
                }
            }

            // Update earliest unmoved index if necessary
            while (earliest_unmoved_index < (A * B) && flip_tracker[earliest_unmoved_index])
                earliest_unmoved_index++;

            // Update bookkeeping for next loop of transpositions
            curr_index = earliest_unmoved_index;
            next_index = get_next_index(curr_index);
            temp_val = data[curr_index];
        }
        return; 
    }
    
    constexpr double TwoPI = 2.0 * M_PI;
    constexpr double NegTwoPI = -TwoPI;
    
    // TODO: Only for circular convolution test! everything 7 seems to all be able to
    //  happen in registers and is sub ns fast!
    void raders_7(float32x2_t *__restrict__ in, float32x2_t *__restrict__ out, float32x2_t *__restrict__ coeff) {
        MCA_START
        //#pragma clang loop unroll_count(7)
        for (unsigned int i = 0; i < 7; i++)
            //#pragma clang loop unroll_count(7)
            for (unsigned int j = 0; j < 7; j++)
                out[i] += neon_mul(in[j], coeff[(j + i) % 7]);
        MCA_END
        return;   
    }

    template <std::size_t N>
    class FFTPlan {
        public:
        static constexpr std::size_t A = prime_factor::get_prime_factor(N);
        static constexpr std::size_t B = N / A;

        FFTPlan() {
            // Ensure the coeffs are calculated
            volatile auto* coefs = gen_coefs_by_angle();
        }

        // TODO : may want to shift the result so that DC component is at the center,
        //  but without this shift the FFT and IFFT are inverses of each other
        static void fft(float32x2_t *__restrict__ data) noexcept {
            fft_recurse(data);

            // Normalize
            for (std::size_t i = 0; i < N; i++)
                data[i] /= static_cast<float>(N);

            return;
        }

        static void ifft(float32x2_t *__restrict__ data) noexcept {
            for (std::size_t i = 0; i < N; i++)
                data[i] = neon_conj(data[i]);

            fft_recurse(data);

            for (std::size_t i = 0; i < N; i++)
                data[i] = neon_conj(data[i]);

            // NOTE : if the fft is normalized, the ifft does not need to be normalized
            //  to keep the ifft the exact inverse operation of the fft.

            return;
        }
        
        // In-place Mixed Radix Cooley Tukey FFT
        static void fft_recurse(float32x2_t *__restrict__ data) noexcept {
            // Base case
            if constexpr (N < 2) {
                return;
            } else if constexpr (N == 2) {
                const float32x2_t a = data[0] + data[1];
                data[1] = data[0] - data[1];
                data[0] = a;
                return;
            }

            std::array<float32x2_t, N> temp;
            if constexpr (B == 1) {
                for (std::size_t i = 0; i < A; i++)
                    temp[i] = data[i];
            } else {
                // Transpose Input Data around the radix
                DFT_binner<B, A>(data, temp.data());

                // Recursively apply fft to each bin
                for (std::size_t i = 0; i < A; i++)
                    FFTPlan<B>::fft_recurse(temp.data() + (i * B));
            }

            // Do tensor multiplication of coeffs with data
            float32x2_t* coefs = gen_coefs_by_angle();
            for (std::size_t i = 0; i < B; i++) {
                #pragma clang loop unroll_count(32)
                for (std::size_t j = 0; j < A; j++) {
                    //MCA_START
                    float32x2_t temp_ans = {0, 0};
                    auto index = i * A * A + j * A;
                    #pragma clang loop unroll_count(32)
                    for (std::size_t k = 0; k < A; k++) {
                        const auto coeff = coefs[index + k];
                        if constexpr (B == 1) {
                            const auto temp_val = temp[k];
                            temp_ans += neon_mul(temp_val, coeff);
                        } else {
                            const auto temp_val = temp[i + (k * B)];
                            temp_ans += neon_mul(temp_val, coeff);
                        }
                    }

                    // Store the answer
                    if constexpr (B == 1) {
                        data[j] = temp_ans;
                    } else {
                        data[i + (j * B)] = temp_ans;
                    }
                    //MCA_END
                }
            }

            // This function is in-place! the input is modified with the result.
            return;
        } 
         
        private:
        // NOTE : coeffs with either method can be calculated as they are needed in the loops that
        //   calculate the actual DFT components, but calculating them ahead of time reduces
        //   done in the DFT itself for the angle based method, and grealy reduces loop
        //   dependencies for the multiplication based method.
        // static float32x2_t* coefs;
        
        // NOTE : Although this is technically correct to calculate coeffs, it is not as
        //   numerically stable as the angle based method.
        static float32x2_t* get_coefs_by_mult() {
            static float32x2_t* coefs = []{
                float32x2_t* result = new float32x2_t[N * A];
            
                // Setup constant factors for incremental rotations by multiplication
                constexpr float small_angle = 
                    static_cast<float>(NegTwoPI / static_cast<double>(N));
                constexpr float big_angle = 
                    static_cast<float>(NegTwoPI / static_cast<double>(A));
                const float32x2_t small_inc = {std::cosf(small_angle), std::sinf(small_angle)};
                const float32x2_t big_inc = {std::cosf(big_angle), std::sinf(big_angle)};
        
                float32x2_t i_factor = {1, 0};
                for (std::size_t i = 0; i < B; i++) {
                    float32x2_t j_factor = i_factor;
                    for (std::size_t j = 0; j < A; j++) {
                        float32x2_t k_factor = {1, 0};
                        for (std::size_t k = 0; k < A; k++) {
                            result[i * A * A + j * A + k] = k_factor;
                            k_factor = neon_mul(k_factor, j_factor);
                        }
                        j_factor = neon_mul(j_factor, big_inc);
                    }
                    i_factor = neon_mul(i_factor, small_inc);
                }
                return result;
            }();
            return coefs;
        }
        
        static float32x2_t gen_factor_by_angle(std::size_t i, std::size_t j, std::size_t k) {
            const double angle = NegTwoPI * static_cast<double>(i * (k + j * B)) / static_cast<double>(N);
            return {static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
        };
        static float32x2_t* gen_coefs_by_angle() {
            static float32x2_t* coefs = []{
                float32x2_t* result = new float32x2_t[N * A];
            
                for (std::size_t i = 0; i < B; i++) {
                    for (std::size_t j = 0; j < A; j++) {
                        for (std::size_t k = 0; k < A; k++) {
                            const auto index = i * A * A + j * A + k;
                            result[index] = gen_factor_by_angle(k, j, i);
                        }
                    }
                }
                return result;
            }();
            return coefs;
        } 
    };

    template <std::size_t N>
    struct FFTPlanBasic {
        static constexpr auto A = prime_factor::get_prime_factor(N);
        static constexpr auto B = N / A;

        // TODO : may want to shift the result so that DC component is at the center,
        //  but without this shift the FFT and IFFT are inverses of each other
        static void fft(std::complex<float> *data) {
            fft_recurse(data);

            // Normalize
            for (std::size_t i = 0; i < N; i++)
                data[i] /= static_cast<float>(N);

            return;
        } 

        // Inverse FFT is just FFT on conj of input, and then conj of output
        static void ifft(std::complex<float> *data) {
            for (std::size_t i = 0; i < N; i++)
                data[i] = std::conj(data[i]);

            fft_recurse(data);

            for (std::size_t i = 0; i < N; i++)
                data[i] = std::conj(data[i]);

            // NOTE : if the fft is normalized, the ifft does not need to be normalized
            //  to keep the ifft the exact inverse operation of the fft.

            return;
        } 

        // In-place Mixed Radix Cooley Tukey FFT
        static void fft_recurse(std::complex<float> *data) {
            // Base case
            if constexpr (A < 2) {
                return;
            }

            // Transpose Input Data around the radix
            std::array<std::complex<float>, N> temp;
            DFT_binner<B, A>(data, temp.data());

            // Recursively apply fft to each bin
            for (std::size_t i = 0; i < A; i++)
                FFTPlanBasic<B>::fft_recurse(temp.data() + (i * B));

            // Combine recursive results
            constexpr auto get_factor = [](std::size_t i, std::size_t j, std::size_t k) -> std::complex<float> {
                const float angle = static_cast<float>(NegTwoPI * static_cast<double>(i * (k + j * B)) / static_cast<double>(N));
                return {std::cosf(angle), std::sinf(angle)};
            };
            
            for (std::size_t i = 0; i < B; i++) {
                // Recursively calculate strided results in temp location
                for (std::size_t j = 0; j < A; j++) {
                    data[i + (j * B)] = {0, 0};
                    for (std::size_t k = 0; k < A; k++) {
                        data[i + (j * B)] += temp[i + (k * B)] * get_factor(k, j, i);
                    }
                }
            }
            // This function is in-place! the input is modified with the result.
            return;
        }
    }; 

    constexpr void wave_gen(std::complex<float> *time_domain, std::complex<float> *freq_domain,
        std::size_t N,
        unsigned int f = 1,
        unsigned int phase = 0,
        unsigned int amp = 1) {
        for (unsigned int i = 0; i < N; i++) {
            const float angle = static_cast<float>(TwoPI) * static_cast<float>((i * f) + phase)
                                  / static_cast<float>(N);
            time_domain[i] += std::complex<float>{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
        }
        const float angle = static_cast<float>(TwoPI) * 
                                  (static_cast<float>(phase) / static_cast<float>(N));
        freq_domain[f] += std::complex<float>{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
    }

    constexpr void wave_gen(float32x2_t *time_domain, float32x2_t *freq_domain,
        std::size_t N,
        unsigned int f = 1,
        unsigned int phase = 0,
        unsigned int amp = 1) {
        for (unsigned int i = 0; i < N; i++) {
            const float angle = static_cast<float>(TwoPI) * static_cast<float>((i * f) + phase)
                                  / static_cast<float>(N);
            time_domain[i] += float32x2_t{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
        }
        const float angle = static_cast<float>(TwoPI) * 
                                  (static_cast<float>(phase) / static_cast<float>(N));
        freq_domain[f] += float32x2_t{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
    }

    constexpr void wave_gen_lcg(std::complex<float> *time_domain, std::complex<float> *freq_domain,
        std::size_t N) {
        for (unsigned int i = 0; i < 13; i++) {
            wave_gen(time_domain, freq_domain, N,
                ((i + 7) * 3) % N,
                ((i + 5) * 11) % N,
                ((i + 11) * 13) % N);
        }
    }

    constexpr void wave_gen_lcg(float32x2_t *time_domain, float32x2_t *freq_domain,
        std::size_t N) {
        for (unsigned int i = 0; i < 13; i++) {
            wave_gen(time_domain, freq_domain, N,
                ((i + 7) * 3) % N,
                ((i + 5) * 11) % N,
                ((i + 11) * 13) % N);
        }
    }
}
