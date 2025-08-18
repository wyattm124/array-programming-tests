#include <array>
#include <complex>
#include <cstddef>
#include <cmath>
#include <arm_neon.h>
#include "prime_factor.hpp"
#include <iostream>
#include <iomanip>

/// TODO: 
/// (1) - Try templated FFT with float32x4_t aligned custom type
/// (2) - implement Bluestiens's for numbers that do not decompse into these factors
/// (3) - implement a top level search for best number to pad an input to for FFT

// For LLVM MCA Analysis
#define MCA_START __asm volatile("# LLVM-MCA-BEGIN");
#define MCA_END __asm volatile("# LLVM-MCA-END");

namespace FFT {
    // These are the factors from highest priority to lowest priority, with
    //  with the highest priority at the lowest index
    constexpr unsigned int factors[] = {8, 4, 6};
    constexpr unsigned int get_best_factor(unsigned int N) {
        for (auto factor : factors) {
            if (N % factor == 0)
                return factor;
        }
        return prime_factor::get_prime_factor(N);
    }

    // Neon intrinsics 
    inline float neon_abs(const float32x2_t &__restrict__ a) noexcept {
        return std::sqrt(a[0]*a[0] + a[1]*a[1]);
    }
    inline float32x2_t neon_conj(const float32x2_t &__restrict__ a) noexcept {
        return {a[0], -a[1]};
    }
    inline float32x2_t neon_mul(const float32x2_t &__restrict__ a, const float32x2_t &__restrict__ b) noexcept {
        return {a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
    }
    inline float32x2_t neon_mul_conj(const float32x2_t &__restrict__ a, const float32x2_t &__restrict__ b_conj) noexcept {
        return {a[0] * b_conj[0] + a[1] * b_conj[1], a[1] * b_conj[0] - a[0] * b_conj[1]};
    } 

    // For DFT Transpositions
    template <unsigned int A, unsigned int B, typename T>
    void DFT_binner(T *__restrict__ in, T *__restrict__ out) noexcept {
        // Base case
        if constexpr (A == 1 || B == 1) {
            std::copy(in, in + A * B, out);
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
        if constexpr (A == 1 || B == 1) {
            return;
        }

        // Have an array to track which elements have already been placed
        std::array<bool, (A * B)> flip_tracker = {false};

        // Given an index n, return the index it should be moved to
        constexpr auto get_next_index = [=](std::size_t n) -> std::size_t {
            const auto i = n % B;
            const auto j = n / B;
            return i * A + j;
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

    template <unsigned int N>
    class FFTPlan {
        public:
        static constexpr unsigned int A = get_best_factor(N);
        static constexpr unsigned int B = N / A;

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
            // Base case + direct DFT cases
            if constexpr (N == 1) {
                // Do nothing, just return
            } else if constexpr (N == 2) {
                const auto temp = data[0] + data[1];
                data[1] = data[0] - data[1];
                data[0] = temp;
            } else if constexpr (N == 3) {
                // TODO : write all math explicitly? See if that gives a speed boost.
                static constexpr float32x2_t third_root = {-0.5f, -0.86660254f};
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                data[0] = x_0 + x_1 + x_2;
                data[1] = x_0 + neon_mul(x_1, third_root) + neon_mul_conj(x_2, third_root);
                data[2] = x_0 + neon_mul_conj(x_1, third_root) + neon_mul(x_2, third_root);
            } else if constexpr (N == 4) {
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                const float32x2_t x_3 = data[3];
                { // First do evens
                    const float32x2_t a = x_0 + x_2;
                    const float32x2_t b = x_1 + x_3;
                    data[0] = a + b;
                    data[2] = a - b;
                }
                { // Then do odds
                    const float32x2_t a = x_0 - x_2;
                    const float32x2_t temp_b = x_1 - x_3;
                    const float32x2_t b = {temp_b[1], -temp_b[0]};
                    data[1] = a + b;
                    data[3] = a - b;
                }
            } else if constexpr (N == 5) {
                static constexpr float32x2_t root_1 = {0.30901699, -0.95105652};
                static constexpr float32x2_t root_2 = {-0.80901699, -0.58778525};
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                const float32x2_t x_3 = data[3];
                const float32x2_t x_4 = data[4];
                data[0] = x_0 + x_1 + x_2 + x_3 + x_4;
                data[1] = x_0 + neon_mul(x_1, root_1) + neon_mul(x_2, root_2) +
                               neon_mul_conj(x_3, root_2) + neon_mul_conj(x_4, root_1);
                data[2] = x_0 + neon_mul(x_1, root_2) + neon_mul_conj(x_2, root_1) +
                               neon_mul(x_3, root_1) + neon_mul_conj(x_4, root_2);
                data[3] = x_0 + neon_mul_conj(x_1, root_2) + neon_mul(x_2, root_1) +
                               neon_mul_conj(x_3, root_1) + neon_mul(x_4, root_2);
                data[4] = x_0 + neon_mul_conj(x_1, root_1) + neon_mul_conj(x_2, root_2) +
                               neon_mul(x_3, root_2) + neon_mul(x_4, root_1);
            } else if constexpr (N == 6) {
                static constexpr float r_a = 0.5f;
                static constexpr float r_b = -0.8660254f;
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                const float32x2_t x_3 = data[3];
                const float32x2_t x_4 = data[4];
                const float32x2_t x_5 = data[5];
                {
                    const float32x2_t a = x_0 + x_2 + x_4;
                    const float32x2_t b = x_1 + x_3 + x_5;
                    data[0] = a + b;
                    data[3] = a - b;
                }
                {
                    const float32x2_t a = x_0 - x_3;
                    const float32x2_t b = x_1 - x_4;
                    const float32x2_t c = x_2 - x_5;
                    data[1] = a + neon_mul(b, float32x2_t{r_a, r_b}) + neon_mul(c, float32x2_t{-r_a, r_b});
                    data[5] = a + neon_mul(b, float32x2_t{r_a, -r_b}) + neon_mul(c, float32x2_t{-r_a, -r_b});
                }
                {
                    const float32x2_t a = x_0 + x_3;
                    const float32x2_t b = x_1 + x_4;
                    const float32x2_t c = x_2 + x_5;
                    data[2] = a + neon_mul(b, float32x2_t{-r_a, r_b}) + neon_mul(c, float32x2_t{-r_a, -r_b});
                    data[4] = a + neon_mul(b, float32x2_t{-r_a, -r_b}) + neon_mul(c, float32x2_t{-r_a, r_b});
                }
            } else if constexpr (N == 7) {
                static constexpr float32x2_t r_1 = {0.62348980, -0.78183148};
                static constexpr float32x2_t r_2 = {-0.22252093, -0.97492791};
                static constexpr float32x2_t r_3 = {-0.90096887, -0.43388374};
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                const float32x2_t x_3 = data[3];
                const float32x2_t x_4 = data[4];
                const float32x2_t x_5 = data[5];
                const float32x2_t x_6 = data[6];
                data[0] = x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6;
                data[1] = x_0 + neon_mul(x_1, r_1) + neon_mul(x_2, r_2) + neon_mul(x_3, r_3) +
                  neon_mul_conj(x_4, r_3) + neon_mul_conj(x_5, r_2) + neon_mul_conj(x_6, r_1);
                data[2] = x_0 + neon_mul(x_1, r_2) + neon_mul_conj(x_2, r_3) + 
                  neon_mul_conj(x_3, r_1) + neon_mul(x_4, r_1) + neon_mul(x_5, r_3) +
                  neon_mul_conj(x_6, r_2);
                data[3] = x_0 + neon_mul(x_1, r_3) + neon_mul_conj(x_2, r_1) + neon_mul(x_3, r_2) +
                  neon_mul_conj(x_4, r_2) + neon_mul(x_5, r_1) + neon_mul_conj(x_6, r_3);
                data[4] = x_0 + neon_mul_conj(x_1, r_3) + neon_mul(x_2, r_1) + 
                  neon_mul_conj(x_3, r_2) + neon_mul(x_4, r_2) + neon_mul_conj(x_5, r_1) +
                  neon_mul(x_6, r_3);
                data[5] = x_0 + neon_mul_conj(x_1, r_2) + neon_mul(x_2, r_3) + neon_mul(x_3, r_1) +
                  neon_mul_conj(x_4, r_1) + neon_mul_conj(x_5, r_3) + neon_mul(x_6, r_2);
                data[6] = x_0 + neon_mul_conj(x_1, r_1) + neon_mul_conj(x_2, r_2) +
                  neon_mul_conj(x_3, r_3) + neon_mul(x_4, r_3) + neon_mul(x_5, r_2) +
                  neon_mul(x_6, r_1);
            } else if constexpr (N == 8) {
                static constexpr float c = 0.70710678118f;
                static constexpr float32x2_t eighth_root = {c,c};
                static constexpr auto forth_root = [](float32x2_t x){
                    return float32x2_t{x[1], -x[0]};
                };
                const float32x2_t x_0 = data[0];
                const float32x2_t x_1 = data[1];
                const float32x2_t x_2 = data[2];
                const float32x2_t x_3 = data[3];
                const float32x2_t x_4 = data[4];
                const float32x2_t x_5 = data[5];
                const float32x2_t x_6 = data[6];
                const float32x2_t x_7 = data[7];
                {
                    const float32x2_t a_0 = x_0 + x_4;
                    const float32x2_t a_1 = x_2 + x_6;
                    const float32x2_t b_0 = x_1 + x_5;
                    const float32x2_t b_1 = x_3 + x_7;
                    { // 0 and 4 index, based on 2 DFTs of size 4 + shifted for index 2 and 6
                        const float32x2_t a = a_0 + a_1;
                        const float32x2_t b = b_0 + b_1;
                        data[0] = a + b;
                        data[4] = a - b;
                    }
                    { // 2 and 6 index, based on 2 DFTs of size 4 + shifted for index 0 and 4
                        const float32x2_t a = a_0 - a_1;
                        const float32x2_t b_p = b_0 - b_1;
                        const float32x2_t b = forth_root(b_p);
                        data[2] = a + b;
                        data[6] = a - b;
                    }
                }
                {
                    const float32x2_t a_0 = x_0 - x_4;
                    const float32x2_t a_1_p = x_2 - x_6;
                    const float32x2_t a_1 = forth_root(a_1_p);
                    { // index 3 and 5
                        const float32x2_t b_0_p = x_1 - x_5;
                        const float32x2_t b_0 = forth_root(b_0_p);
                        const float32x2_t b_1 = x_3 - x_7;
                        data[3] = (a_0 - a_1) + neon_mul_conj(b_1 + b_0, eighth_root);
                        data[5] = (a_0 + a_1) + neon_mul(b_1 - b_0, eighth_root);
                    }
                    { // DFT for index 1 and 7
                        const float32x2_t b_0 = x_1 - x_5;
                        const float32x2_t b_1_p = x_3 - x_7;
                        const float32x2_t b_1 = forth_root(b_1_p);
                        data[1] = (a_0 + a_1) + neon_mul_conj(b_0 + b_1, eighth_root);
                        data[7] = (a_0 - a_1) + neon_mul(b_0 - b_1, eighth_root); 
                    }
                }
            } else if constexpr (A == 1 || B == 1) { // Do the full DFT
                auto dft_matrix = gen_coefs_by_angle();
                std::array<float32x2_t, N> temp_data;

                // Do a full DFT
                for (unsigned int i = 0; i < N; i++) {
                    float32x2_t temp_ans = {0, 0};
                    for (unsigned int j = 0; j < N; j++) {
                        temp_ans += neon_mul(data[j], dft_matrix[i * N + j]);
                    }
                    temp_data[i] = temp_ans;
                }

                for (unsigned int i = 0; i < N; i++)
                    data[i] = temp_data[i];
            } else { // Do a mixed radix FFT
                std::array<float32x2_t, N> temp_data;

                // Transpose Input Data around the first radix
                DFT_binner<A, B>(data, temp_data.data());
 
                // Do the A sized FFT on each bin
                for (unsigned int i = 0; i < B; i++)
                    FFTPlan<A>::fft_recurse(temp_data.data() + (i * A));

                // Multiply by Cooley Tukey Twiddle factors
                auto twiddle_factors = gen_coefs_by_angle();
                for (unsigned int i = 0; i < N; i++) {
                    data[i] = neon_mul(temp_data[i], twiddle_factors[((i / A) * (i % A)) % N]);
                }

                // Transpose around the second radix
                DFT_binner<B, A>(data, temp_data.data());

                // Do the B sized FFT on each bin
                for (unsigned int i = 0; i < A; i++)
                    FFTPlan<B>::fft_recurse(temp_data.data() + (i * B));

                DFT_binner<A, B>(temp_data.data(), data);
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
              if constexpr (A == 1 || B == 1) { // DFT matrix
                float32x2_t* result = new float32x2_t[N * N];
                for (std::size_t i = 0; i < N; i++) {
                  for (std::size_t j = 0; j < N; j++) {
                    const double angle = 
                      NegTwoPI * static_cast<double>((i * j) % N) / static_cast<double>(N);
                        result[i * N + j] = float32x2_t{static_cast<float>(std::cos(angle)),
                                                        static_cast<float>(std::sin(angle))};
                  }
                }
                return result;
              } else { // Split Radix FFT Twiddle factors
                float32x2_t* result = new float32x2_t[N];
                for (std::size_t i = 0; i < N; i++) {
                  const double angle = 
                      NegTwoPI * static_cast<double>(i) / static_cast<double>(N);
                  result[i] = {static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
                }
                return result;
              }
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

    template <typename T>
    constexpr void wave_gen(T *time_domain, T *freq_domain,
        std::size_t N,
        unsigned int f = 1,
        unsigned int phase = 0,
        unsigned int amp = 1) {
        for (unsigned int i = 0; i < N; i++) {
            const float angle = static_cast<float>(TwoPI) * static_cast<float>((i * f) + phase)
                                  / static_cast<float>(N);
            time_domain[i] += T{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
        }
        const float angle = static_cast<float>(TwoPI) * 
                                  (static_cast<float>(phase) / static_cast<float>(N));
        freq_domain[f] += T{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
    } 

    template <typename T>
    constexpr void wave_gen_lcg(T *time_domain, T *freq_domain, unsigned int N) {
        if (N < 13) {
            if (N > 7) {
                wave_gen(time_domain, freq_domain, N,
                    7, 2, 1);
            }
            if (N > 5) {
                wave_gen(time_domain, freq_domain, N,
                    5, 2, 1);
            }
            if (N > 4) {
                wave_gen(time_domain, freq_domain, N,
                    4, 3, 2);
            }
            if (N > 3) {
                wave_gen(time_domain, freq_domain, N,
                    3, 2, 1);
            }
            wave_gen(time_domain, freq_domain, N,
                1, 1, 1);
            return;
        } else {
            for (unsigned int i = 0; i < 13; i++) {
                wave_gen(time_domain, freq_domain, N,
                    ((i + 7) * 3) % N,
                    ((i + 5) * 11) % N,
                    ((i + 11) * 13) % N);
            }
        }
    }
}
