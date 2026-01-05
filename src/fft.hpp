#pragma once

#include <cmath>
#include <array>
#include <new>
#include <memory>

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

// For Cache aligning
#define MY_CPU_LOAD_SIZE 64

namespace FFT {

    /// Signal generation for generating an FFT input with the corresponding expected output.
    constexpr double TwoPI = 2.0 * M_PI;
    
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

    /// FFT factor generation
    template<typename T>
    struct FFTFactorGen {
        // For this type parameter, these are the factors, from highest priority to
        //  lowest priority, to create each layer of the FFT.
        static constexpr unsigned int factors[] = {8, 4, 6};

        // Given an input of size N, return the FFT bin sizes of the next
        //  layer of the FFT.
        static constexpr unsigned int get_best_factor(unsigned int N) {
            for (auto factor : factors) {
                if (N % factor == 0)
                    return factor;
            }
            return prime_factor::get_prime_factor(N);
        }
    };

    /// FFT coefficient generation
    constexpr double NegTwoPI = -TwoPI;

    template<typename T, unsigned int N>
    void populate_dft_matrix_by_angle(T *mat) {
        for (std::size_t i = 0; i < N; i++) {
            for (std::size_t j = 0; j < N; j++) {
                const double angle = 
                    NegTwoPI * static_cast<double>((i * j) % N) / static_cast<double>(N);
                    mat[i * N + j] = T{static_cast<float>(std::cos(angle)),
                                          static_cast<float>(std::sin(angle))};
            }
        }
    }
    
    template<typename T, unsigned int A, unsigned int B>
    void populate_twiddle_factors_by_angle(T* factors) {
        for (unsigned int i = 0; i < A; i++) {
            for (unsigned int j = 0; j < B; j++) {
                const double angle = NegTwoPI * static_cast<double>(i * j) / static_cast<double>(A * B);
                factors[i * B + j] = T{static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
            }
        }
    } 
    
    template <typename T>
    class FFTPlan {
        public:
        // Multiplication op aliases for type parameter
        static constexpr auto m = mult<T>;
        static constexpr auto mc = mult_conj<T>;

        template<unsigned int N>
        static void Init() { FFTLayer<N>::Init(); }

        // The result is not shifted so the DC component is still at index 0,
        //  and the FFT and IFFT are inverses of each other
        template<unsigned int N>
        static void fft(T *__restrict__ in, T *__restrict__ out) noexcept {
            if constexpr (FFTLayer<N>::base_case) {
                FFTLayer<N>::fft_recurse(in, out);
            } else {
                alignas(MY_CPU_LOAD_SIZE) std::array<T, N> temp_data;
                FFTLayer<N>::fft_recurse(in, temp_data.data());
                transpose<N, 1>(temp_data.data(), out);
            }

            // Normalize
            // FFTW does not normalize, so this loop should be commented
            //  out for appropriate performance comparisons.
            for (std::size_t i = 0; i < N; i++)
                out[i] /= static_cast<float>(N);

            return;
        }

        // The Inverse DFT matrix is the same as the DFT matrix but with the corresponding
        //  factors conjugated.
        template<unsigned int N>
        static void ifft(T *__restrict__ in, T *__restrict__ out) noexcept {
            for (std::size_t i = 0; i < N; i++)
                in[i] = conj(in[i]);
            
            if constexpr (FFTLayer<N>::base_case) {
                FFTLayer<N>::fft_recurse(in, out);
            } else {
                alignas(MY_CPU_LOAD_SIZE) std::array<T, N> temp_data;
                FFTLayer<N>::fft_recurse(in, temp_data.data());
                transpose<N, 1>(temp_data.data(), out);
            }

            for (std::size_t i = 0; i < N; i++)
                out[i] = conj(out[i]);

            // NOTE : if the fft is normalized, the ifft does not need to be normalized
            //  to keep the ifft the exact inverse operation of the fft.

            return;
        }

        private:

        template<unsigned int N, unsigned int S>
        static void transpose(T *__restrict__ in, T *__restrict__ out) noexcept {
            constexpr unsigned int A = FFTFactorGen<T>::get_best_factor(N);
            constexpr unsigned int B = N / A;

            // Transpose result back in order
            for (unsigned int i = 0; i < A; i++) {
                if constexpr (B == 1) {
                    out[S * i] = in[i];
                } else {
                    transpose<B, A * S>(in + (B * i), out + (S * i));
                }
            }
            return;
        }

        template<unsigned int N>
        struct FFTLayer {
            static constexpr unsigned int A = FFTFactorGen<T>::get_best_factor(N);
            static constexpr unsigned int B = N / A;
            static constexpr bool base_case = false;

            static void Init() {
                // Ensure the coeffs are calculated
                volatile T* coefs = FFTPlan<T>::get_twiddle_factors_by_angle<A, B>();
                FFTLayer<A>::Init();
                FFTLayer<B>::Init();

                // First layer for A must be a non-recursive base case to make sure
                //  recursive transpositions are done as expected
                static_assert(FFTLayer<A>::base_case);
            }

            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                // TODO : may not want to put this on the stack for large FFT steps
                alignas(MY_CPU_LOAD_SIZE) std::array<T, N> temp_data;

                for (unsigned int i = 0; i < B; i++) {

                    alignas(MY_CPU_LOAD_SIZE) std::array<T, A> workspace;

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
                    FFTLayer<A>::fft_recurse(workspace.data(), temp_data.data() + (i * A));
                }

                // Cooley Tukey twiddle factors
                T *twiddle_factors = std::assume_aligned<MY_CPU_LOAD_SIZE>(
                    FFTPlan<T>::get_twiddle_factors_by_angle<A, B>());

                for (unsigned int i = 0; i < A; i++) {

                    alignas(MY_CPU_LOAD_SIZE) std::array<T, B> workspace;
                 
                    /* Same tactic as above to transpose input data around
                     *  second radix B.
                     * 
                     * Here is different than above only in tha before using
                     *  this input for the recursive FFT, we make sure to
                     *  adjust it by the appropriate twiddle factor.
                     */

                    // TODO : Need to order twiddle factors so they are accessed linearly
                    for (unsigned int j = 0; j < B; j++) {
                        workspace[j] = m(temp_data[i + j * A], twiddle_factors[i * B + j]);
                    }

                    // Do the B sized FFT on the current bin
                    FFTLayer<B>::fft_recurse(workspace.data(), out + (i * B));
                }
            }
        };

        // Partial specialization for prime N - uses full DFT matrix
        template<unsigned int N>
        requires prime_factor::Prime<N>
        struct FFTLayer<N> {
            static constexpr unsigned int A = N;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {
                // Ensure the coeffs are calculated
                volatile T* coefs = FFTPlan<T>::get_dft_matrix_by_angle<N>();
            }
            
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                static T* dft_matrix = std::assume_aligned<MY_CPU_LOAD_SIZE>(
                    FFTPlan<T>::get_dft_matrix_by_angle<N>());

                for (unsigned int i = 0; i < N; i++) {
                    T temp_ans = {0, 0};
                    for (unsigned int j = 0; j < N; j++) {
                        temp_ans += m(in[j], dft_matrix[i * N + j]);
                    }
                    out[i] = temp_ans;
                }
            }
        };

        // Hand written DFT cases serve as base cases for recursive Cooley Tookey
        template<>
        struct FFTLayer<1> {
            static constexpr unsigned int A = 1;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {
                // This init should never be called, as this FFTLayer should never be used
                static_assert(false);
            }
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                // All base cases should avoid the need for this trivial specialization
                static_assert(false);
            }
        };
        
        template<>
        struct FFTLayer<2> {
            static constexpr unsigned int A = 2;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                const auto x_0 = in[0];
                const auto x_1 = in[1];
                out[0] = x_0 + x_1;
                out[1] = x_0 - x_1;
                return;
            }
        };

        template<>
        struct FFTLayer<3> {
            static constexpr unsigned int A = 3;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr T third_root = {-0.5f, -0.866025403784f};
                const T x_0 = in[0];
                const T x_1 = in[1];
                const T x_2 = in[2];
                out[0] = x_0 + x_1 + x_2;
                out[1] = x_0 + m(x_1, third_root) + mc(x_2, third_root);
                out[2] = x_0 + mc(x_1, third_root) + m(x_2, third_root);
                return;
            }
        };

        template<>
        struct FFTLayer<4> {
            static constexpr unsigned int A = 4;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
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
                return;
            }
        };

        template<>
        struct FFTLayer<5> {
            static constexpr unsigned int A = 5;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
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
                return;
            }
        };

        template<>
        struct FFTLayer<6> {
            static constexpr unsigned int A = 6;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
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
                return;
            }
        };

        template<>
        struct FFTLayer<7> {
            static constexpr unsigned int A = 7;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
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
                return;
            }
        };

        template<>
        struct FFTLayer<8> {
            static constexpr unsigned int A = 8;
            static constexpr unsigned int B = 1;
            static constexpr bool base_case = true;
            static void Init() {}
            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
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
                return;
            }
        };

        private:
        // Coeffs calculated ahead of time trades memory for speed.
        //   Since Coeffs are roots of unity they could also be multiplied
        //   together within loops to generate them as needed, but this 
        //   has emperically been found to be not acceptably stable numerically.
        //   It is better to generate them with this angle based method for 
        //   numerical stablility. 

        // TODO : use a cleaner way to hold coefs at class rather than object level
        template<unsigned int N>
        static T* get_dft_matrix_by_angle() {
            static T* coefs = []{
                T* result = static_cast<T*>(
                    ::operator new(sizeof(T) * N * N, std::align_val_t(MY_CPU_LOAD_SIZE)));
                populate_dft_matrix_by_angle<T, N>(result); 
                return result;
            }();
            return coefs;
        };

        template<unsigned int A, unsigned int B>
        static T* get_twiddle_factors_by_angle() {
            static T* coefs = []{
                T* result = static_cast<T*>(
                    ::operator new(sizeof(T) * A * B, std::align_val_t(MY_CPU_LOAD_SIZE)));
                populate_twiddle_factors_by_angle<T, A, B>(result); 
                return result;
            }();
            return coefs;
        };
    }; 
}
