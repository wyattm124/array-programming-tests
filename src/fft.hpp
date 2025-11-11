#pragma once

#include <cmath>
#include <array>

#include "prime_factor.hpp"

/// TODO:
/// (1) - May want base cases for 9 and 10
/// (2) - Implement Rader's algorithm for small prime numbers, especially 11, 13, 17 and 19
/// (3) - Implement Bluestiens's for numbers that do not decompse into any hand written cases
/// (4) - Implement a top level search to appropriately pad an input to a number if needed
///   and factor it into the best radices

// For LLVM MCA Analysis
#define MCA_START __asm volatile("# LLVM-MCA-BEGIN");
#define MCA_END __asm volatile("# LLVM-MCA-END");

namespace FFT {
    /// Operations used by FFT
    template <typename T>
    inline float abs(const T &__restrict__ a) noexcept {
        return std::sqrt(a[0]*a[0] + a[1]*a[1]);
    }

    template <typename T>
    constexpr T conj(const T &__restrict__ a) noexcept {
        return {a[0], -a[1]};
    }

    template <typename T>
    constexpr T flipper(const T &__restrict__ a) noexcept {
        return {a[1], a[0]};
    }

    template <typename T>
    inline T mult(const T &__restrict__ a, const T &__restrict__ b) noexcept {
        return {a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
    }
    template <typename T>
    inline T mult_conj(const T &__restrict__ a, const T &__restrict__ b_conj) noexcept {
        return {a[0] * b_conj[0] + a[1] * b_conj[1], a[1] * b_conj[0] - a[0] * b_conj[1]};
    }

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

    // For DFT Transpositions
    template <unsigned int A, unsigned int B, typename T>
    void DFT_binner(T *__restrict__ in, T *__restrict__ out) noexcept {
        // Base case
        if constexpr (A == 1 || B == 1) {
            for (unsigned int i = 0; i < A * B; i++)
                out[i] = in[i];
            return;
        }

        // This will create B DFTs of size A
        for (unsigned int i = 0; i < B; i++) {
            for (unsigned int j = 0; j < A; j++) {
                out[i * A + j] = in[j * B + i];
            }
        }
        return;
    }
    
    // NOTE : although this algorithm is in-place, all the added computation is not worth
    //  the memory savings.

    // This will create B DFTs of size A
    template <unsigned int A, unsigned int B, typename T>
    void DFT_binner_inplace(T *__restrict__ data) {
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
        T temp_val = data[curr_index];
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

    template <unsigned int N, typename T>
    class FFTPlan {
        private:
        static constexpr auto m = mult<T>;
        static constexpr auto mc = mult_conj<T>;
        public:
        static constexpr unsigned int A = get_best_factor(N);
        static constexpr unsigned int B = N / A; 

        FFTPlan() {
            // Ensure the coeffs are calculated
            volatile auto* coefs = gen_coefs_by_angle();
        }

        // NOTE : The result is not shifted so the DC component is still at index 0,
        //  and the FFT and IFFT are inverses of each other
        static void fft(T *__restrict__ in, T *__restrict__ out) noexcept {
            fft_recurse(in, out);

            // Normalize
            for (std::size_t i = 0; i < N; i++)
                out[i] /= static_cast<float>(N);

            return;
        }

        // The Inverse DFT matrix is the same as the DFT matrix but with the corresponding
        //  factors conjugated.
        static void ifft(T *__restrict__ in, T *__restrict__ out) noexcept {
            for (std::size_t i = 0; i < N; i++)
                in[i] = conj(in[i]);

            fft_recurse(in, out);

            for (std::size_t i = 0; i < N; i++)
                out[i] = conj(out[i]);

            // NOTE : if the fft is normalized, the ifft does not need to be normalized
            //  to keep the ifft the exact inverse operation of the fft.

            return;
        }
        
        // Mixed Radix Cooley Tukey FFT
        static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
            MCA_START
            // Base case + direct DFT cases
            if constexpr (N == 1) {
                out[0] = in[0];
            } else if constexpr (N == 2) {
                const auto x_0 = in[0];
                const auto x_1 = in[1];
                out[0] = x_0 + x_1;
                out[1] = x_0 - x_1;
            } else if constexpr (N == 3) {
                // TODO : write all math explicitly? See if that gives a speed boost.
                static constexpr T third_root = {-0.5f, -0.866025403784f};
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                out[0] = x_0 + x_1 + x_2;
                out[1] = x_0 + m(x_1, third_root) + mc(x_2, third_root);
                out[2] = x_0 + mc(x_1, third_root) + m(x_2, third_root);
            } else if constexpr (N == 4) {
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                const T x_3 = in[3];
                { // First do evens
                    const T a = x_0 + x_2;
                    const T b = x_1 + x_3;
                    out[0] = a + b;
                    out[2] = a - b;
                }
                { // Then do odds
                    const T a = x_0 - x_2;
                    const T temp_b = x_1 - x_3;
                    const T b = {temp_b[1], -temp_b[0]};
                    out[1] = a + b;
                    out[3] = a - b;
                }
            } else if constexpr (N == 5) {
                static constexpr T root_1 = {0.309016994375, -0.951056516295};
                static constexpr T root_2 = {-0.809016994375, -0.587785252292};
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                const T x_3 = in[3];
                const T x_4 = in[4];
                out[0] = x_0 + x_1 + x_2 + x_3 + x_4;
                out[1] = x_0 + m(x_1, root_1) + m(x_2, root_2) +
                               mc(x_3, root_2) + mc(x_4, root_1);
                out[2] = x_0 + m(x_1, root_2) + mc(x_2, root_1) +
                               m(x_3, root_1) + mc(x_4, root_2);
                out[3] = x_0 + mc(x_1, root_2) + m(x_2, root_1) +
                               mc(x_3, root_1) + m(x_4, root_2);
                out[4] = x_0 + mc(x_1, root_1) + mc(x_2, root_2) +
                               m(x_3, root_2) + m(x_4, root_1);
            } else if constexpr (N == 6) {
                static constexpr float r_a = 0.5f;
                static constexpr float r_b = -0.866025403784f;
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                const T x_3 = in[3];
                const T x_4 = in[4];
                const T x_5 = in[5];
                {
                    const T a = x_0 + x_2 + x_4;
                    const T b = x_1 + x_3 + x_5;
                    out[0] = a + b;
                    out[3] = a - b;
                }
                {
                    const T a = x_0 - x_3;
                    const T b = x_1 - x_4;
                    const T c = x_2 - x_5;
                    out[1] = a + m(b, T{r_a, r_b}) + m(c, T{-r_a, r_b});
                    out[5] = a + m(b, T{r_a, -r_b}) + m(c, T{-r_a, -r_b});
                }
                {
                    const T a = x_0 + x_3;
                    const T b = x_1 + x_4;
                    const T c = x_2 + x_5;
                    out[2] = a + m(b, T{-r_a, r_b}) + m(c, T{-r_a, -r_b});
                    out[4] = a + m(b, T{-r_a, -r_b}) + m(c, T{-r_a, r_b});
                }
            } else if constexpr (N == 7) {
                static constexpr T r_1 = {0.623489801859, -0.781831482468};
                static constexpr T r_2 = {-0.222520933956, -0.974927912182};
                static constexpr T r_3 = {-0.900968867902, -0.433883739118};
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                const T x_3 = in[3];
                const T x_4 = in[4];
                const T x_5 = in[5];
                const T x_6 = in[6];
                out[0] = x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6;
                out[1] = x_0 + m(x_1, r_1) + m(x_2, r_2) + m(x_3, r_3) +
                  mc(x_4, r_3) + mc(x_5, r_2) + mc(x_6, r_1);
                out[2] = x_0 + m(x_1, r_2) + mc(x_2, r_3) + 
                  mc(x_3, r_1) + m(x_4, r_1) + m(x_5, r_3) +
                  mc(x_6, r_2);
                out[3] = x_0 + m(x_1, r_3) + mc(x_2, r_1) + m(x_3, r_2) +
                  mc(x_4, r_2) + m(x_5, r_1) + mc(x_6, r_3); 
                out[4] = x_0 + mc(x_1, r_3) + m(x_2, r_1) + 
                  mc(x_3, r_2) + m(x_4, r_2) + mc(x_5, r_1) +
                  m(x_6, r_3);
                out[5] = x_0 + mc(x_1, r_2) + m(x_2, r_3) + m(x_3, r_1) +
                  mc(x_4, r_1) + mc(x_5, r_3) + m(x_6, r_2);
                out[6] = x_0 + mc(x_1, r_1) + mc(x_2, r_2) +
                  mc(x_3, r_3) + m(x_4, r_3) + m(x_5, r_2) +
                  m(x_6, r_1);
            } else if constexpr (N == 8) {
                static constexpr float c = 0.707106781187f;
                static constexpr T eighth_root = {c,c};
                static constexpr auto forth_root = [](T x){
                    return T{x[1], -x[0]};
                };
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                const T x_3 = in[3];
                const T x_4 = in[4];
                const T x_5 = in[5];
                const T x_6 = in[6];
                const T x_7 = in[7];
                {
                    const T a_0 = x_0 + x_4;
                    const T a_1 = x_2 + x_6;
                    const T b_0 = x_1 + x_5;
                    const T b_1 = x_3 + x_7;
                    { // 0 and 4 index, based on 2 DFTs of size 4 + shifted for index 2 and 6
                        const T a = a_0 + a_1;
                        const T b = b_0 + b_1;
                        out[0] = a + b;
                        out[4] = a - b;
                    }
                    { // 2 and 6 index, based on 2 DFTs of size 4 + shifted for index 0 and 4
                        const T a = a_0 - a_1;
                        const T b_p = b_0 - b_1;
                        const T b = forth_root(b_p);
                        out[2] = a + b;
                        out[6] = a - b;
                    }
                }
                {
                    const T a_0 = x_0 - x_4;
                    const T a_1_p = x_2 - x_6;
                    const T a_1 = forth_root(a_1_p);
                    { // index 3 and 5
                        const T b_0_p = x_1 - x_5;
                        const T b_0 = forth_root(b_0_p);
                        const T b_1 = x_3 - x_7;
                        out[3] = (a_0 - a_1) + mc(b_1 + b_0, eighth_root);
                        out[5] = (a_0 + a_1) + m(b_1 - b_0, eighth_root);
                    }
                    { // DFT for index 1 and 7
                        const T b_0 = x_1 - x_5;
                        const T b_1_p = x_3 - x_7;
                        const T b_1 = forth_root(b_1_p);
                        out[1] = (a_0 + a_1) + mc(b_0 + b_1, eighth_root);
                        out[7] = (a_0 - a_1) + m(b_0 - b_1, eighth_root); 
                    }
                }
            } else if constexpr (A == 1 || B == 1) { // Do the full DFT
                auto dft_matrix = gen_coefs_by_angle();
                std::array<T, N> temp_data;

                for (unsigned int i = 0; i < N; i++) {
                    T temp_ans = {0, 0};
                    for (unsigned int j = 0; j < N; j++) {
                        temp_ans += m(in[j], dft_matrix[i * N + j]);
                    }
                    out[i] = temp_ans;
                }
            } else { // Do a recursive step of a mixed radix FFT
                // TODO : may not want to put this on the stack for large FFT steps
                alignas(16) std::array<T, N> temp_data;

                for (unsigned int i = 0; i < B; i++) {
                    alignas(16) std::array<T, A> workspace;

                    /* Transpose input data around the first radix A
                     *  to create B bins of size A.
                     * 
                     * We create a single bin and compute its FFT before
                     *  moving onto the next bin to improve cache locality.
                     */
                    for (unsigned int j = 0; j < A; j++) {
                        workspace[j] = in[i + j * B];
                    }

                    // Do the A sized FFT on the current bin
                    FFTPlan<A, T>::fft_recurse(workspace.data(), out + (i * A));
                }

                for (unsigned int i = 0; i < A; i++) {
                    // Cooley Tukey twiddle factors
                    T *twiddle_factors = gen_coefs_by_angle();

                    alignas(16) std::array<T, B> workspace;
                 
                    /* Same tactic as above to transpose input data around
                     *  second radix B.
                     * 
                     * Here is different than above only in tha before using
                     *  this input for the recursive FFT, we make sure to
                     *  adjust it by the appropriate twiddle factor.
                     */
                    for (unsigned int j = 0; j < B; j++) {
                        workspace[j] = m(out[i + j * A], twiddle_factors[i * j]);
                    }

                    // Do the B sized FFT on the current bin
                    FFTPlan<B, T>::fft_recurse(workspace.data(), temp_data.data() + (i * B));
                }

                // Transpose result back in order
                DFT_binner<A, B>(temp_data.data(), out);
            }
            MCA_END

            return;
        } 
         
        private:
        // NOTE : coeffs with either method can be calculated as they are needed in the loops that
        //   calculate the actual DFT components, but calculating them ahead of time reduces
        //   work done in the DFT itself for the angle based method, and grealy reduces loop
        //   dependencies for the multiplication based method.
        // static float32x2_t* coefs;
        
        // NOTE : Although this is technically correct to calculate coeffs, it is not as
        //   numerically stable as the angle based method.
        static T* get_coefs_by_mult() {
            static T* coefs = []{
                T* result = new T[N * A];
            
                // Setup constant factors for incremental rotations by multiplication
                constexpr float small_angle = 
                    static_cast<float>(NegTwoPI / static_cast<double>(N));
                constexpr float big_angle = 
                    static_cast<float>(NegTwoPI / static_cast<double>(A));
                const T small_inc = {std::cosf(small_angle), std::sinf(small_angle)};
                const T big_inc = {std::cosf(big_angle), std::sinf(big_angle)};
        
                T i_factor = {1, 0};
                for (std::size_t i = 0; i < B; i++) {
                    T j_factor = i_factor;
                    for (std::size_t j = 0; j < A; j++) {
                        T k_factor = {1, 0};
                        for (std::size_t k = 0; k < A; k++) {
                            result[i * A * A + j * A + k] = k_factor;
                            k_factor = m(k_factor, j_factor);
                        }
                        j_factor = m(j_factor, big_inc);
                    }
                    i_factor = m(i_factor, small_inc);
                }
                return result;
            }();
            return coefs;
        }
        
        static T gen_factor_by_angle(std::size_t i, std::size_t j, std::size_t k) {
            const double angle = NegTwoPI * static_cast<double>(i * (k + j * B)) / static_cast<double>(N);
            return {static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
        }

        // TODO : use a cleaner way to hold coefs at class rather than object level
        static T* gen_coefs_by_angle() {
            static T* coefs = []{
              if constexpr (A == 1 || B == 1) { // DFT matrix
                T* result = new T[N * N];
                for (std::size_t i = 0; i < N; i++) {
                  for (std::size_t j = 0; j < N; j++) {
                    const double angle = 
                      NegTwoPI * static_cast<double>((i * j) % N) / static_cast<double>(N);
                        result[i * N + j] = T{static_cast<float>(std::cos(angle)),
                                              static_cast<float>(std::sin(angle))};
                  }
                }
                return result;
              } else { // Split Radix FFT Twiddle factors
                T* result = new T[N];
                for (std::size_t i = 0; i < N; i++) {
                  const double angle = 
                      NegTwoPI * static_cast<double>(i) / static_cast<double>(N);
                  result[i] = T{static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
                }
                return result;
              }
            }();
            return coefs;
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
            if (N > 7) wave_gen(time_domain, freq_domain, N, 7, 2, 1);
            if (N > 5) wave_gen(time_domain, freq_domain, N, 5, 2, 1);
            if (N > 4) wave_gen(time_domain, freq_domain, N, 4, 3, 2);
            if (N > 3) wave_gen(time_domain, freq_domain, N, 3, 2, 1);
            wave_gen(time_domain, freq_domain, N, 1, 1, 1);
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
