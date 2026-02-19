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
        for (unsigned int j = 0; j < B; j++) {
            for (unsigned int i = 0; i < A; i++) {
                const double angle = NegTwoPI * static_cast<double>(i * j) / static_cast<double>(A * B);
                factors[j * A + i] = T{static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
            }
        }
    }
    
    template <typename T>
    class FFTPlan {
        public:
        // Multiplication op aliases for type parameter
        static constexpr auto m = mult<T>;
        static constexpr auto mc = mult_conj<T>;

        // Layers need to be initialized recursively to populate
        //  their corresponding coefficient arrays.
        template<unsigned int N>
        static void Init() { FFTRecurseLayer<N,N>::Init(); }

        template<unsigned int N>
        static void fft(T *__restrict__ in, T *__restrict__ out) noexcept {
            top_level_fft<N, true>(in, out);
            return;
        }
 
        template<unsigned int N>
        static void ifft(T *__restrict__ in, T *__restrict__ out) noexcept {
            top_level_fft<N, false>(in, out);
            return;
        }

        private:

        // switchable forward / backward fft top level function to reduce redundant code
        template<unsigned int N, bool forward>
        static void top_level_fft(T *__restrict__ in, T *__restrict__ out) noexcept {
            // The result is not shifted so the DC component is still at index 0,
            //  and the FFT and IFFT are inverses of each other
            
            // The Inverse DFT matrix is the same as the DFT matrix but with the corresponding
            //  factors conjugated.

            if constexpr (FFTRecurseLayer<N,N>::base_case) {
                // If we are already at a base case, then we just need to do a straightforward
                //  DFT.
                DFT<N, forward>::execute(in, out);

                for (int i = 0; i < N; i++) {
                    out[i] = forward ? (out[i] / static_cast<float>(N)) : conj(out[i]);
                }
            } else {
                // TODO : large FFTs overflow the stack, this array may need to be put on the heap.
                alignas(MY_CPU_LOAD_SIZE) std::array<T, N> temp;
                constexpr unsigned int A = FFTRecurseLayer<N,N>::A;
                constexpr unsigned int B = FFTRecurseLayer<N,N>::B;
                for (unsigned int i = 0; i < B; i++) {
                    /* Transpose input data around the first radix A
                     *  to create B bins of size A.
                     *
                     * Generally, in the recursive case each layer is
                     *  execute as:
                     *  (1) Transpose around Radix A
                     *  (2) Compute A sized DFT on each bin
                     *  (3) Transpose around Radix B
                     *  (4) Compute B sized FFT on each bin
                     *    - The NEXT layer's N is the CURRENT layer's B
                     *    - If the NEXT layer's B is 1, this is a base case
                     *      and the layer will just compute the B sized DFT
                     *      on each bin.
                     *    - If the NEXT layer's B is NOT 1, then it will
                     *      be another recursive layer with the new
                     *      N, A, and B.
                     *
                     * The recursive nature of the transpositions allows them to
                     *  be fused between layers. However, this fusion also requires
                     *  this required the first transposition to be completed here.
                     *
                     * To minimize data movement, this initial transposition is also
                     *  used to load the input into an array the FFT can write over
                     *  as workspace memory.
                     */
                        
                    for (unsigned int j = 0; j < A; j++) {
                        out[i * A + j] = forward ? in[i + j * B] : conj(in[i + j * B]);
                    } 
                }
                FFTRecurseLayer<N,N>::fft_recurse(out, temp.data());
                transpose<N, 1, forward>(temp.data(), out);
            } 

            return;
        }
        
        // A full DFT done by simple matrix multiplication with the appropriate DFT matrix.
        //  Specific sizes are specialized and optimized as appropriate.
        template<unsigned int N, bool forward> 
        struct DFT {
            static void Init() {
                // Ensure the coeffs are calculated
                volatile T* coefs = FFTPlan<T>::get_dft_matrix_by_angle<N>();
            }
            
            // Simply compute the multiplication of the input array (vector) by
            //  by the DFT matrix.
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static T* dft_matrix = std::assume_aligned<MY_CPU_LOAD_SIZE>(
                    FFTPlan<T>::get_dft_matrix_by_angle<N>());

                for (unsigned int i = 0; i < N; i++) {
                    T temp_ans = {0, 0};
                    for (unsigned int j = 0; j < N; j++) {
                        temp_ans += m(forward ? in[j] : conj(in[j]), dft_matrix[i * N + j]);
                    }
                    out[i] = temp_ans;
                }
            }
        };

        // Hand written DFT cases serve as base cases for recursive Cooley Tookey
        template<bool forward>
        struct DFT<1, forward> {
            static void Init() {
                // This init should never be called, as this FFTRecurseLayer should never be used
                static_assert(false);
            }
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                // All base cases should avoid the need for this trivial specialization
                static_assert(false);
            }
        };
        
        template<bool forward>
        struct DFT<2, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                const auto x_0 = forward ? in[0] : conj(in[0]);
                const auto x_1 = forward ? in[1] : conj(in[1]);
                out[0] = x_0 + x_1;
                out[1] = x_0 - x_1;
                return;
            }
        };

        template<bool forward>
        struct DFT<3, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr T third_root = {-0.5f, -0.866025403784f};
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                out[0] = x_0 + x_1 + x_2;
                out[1] = x_0 + m(x_1, third_root) + mc(x_2, third_root);
                out[2] = x_0 + mc(x_1, third_root) + m(x_2, third_root);
                return;
            }
        };

        template<bool forward>
        struct DFT<4, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                const T x_3 = forward ? in[3] : conj(in[3]);
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

        template<bool forward>
        struct DFT<5, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr T root_1 = {0.309016994375, -0.951056516295};
                static constexpr T root_2 = {-0.809016994375, -0.587785252292};
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                const T x_3 = forward ? in[3] : conj(in[3]);
                const T x_4 = forward ? in[4] : conj(in[4]);
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

        template<bool forward>
        struct DFT<6, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr float r_a = 0.5f;
                static constexpr float r_b = -0.866025403784f;
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                const T x_3 = forward ? in[3] : conj(in[3]);
                const T x_4 = forward ? in[4] : conj(in[4]);
                const T x_5 = forward ? in[5] : conj(in[5]);
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

        template<bool forward>
        struct DFT<7, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr T r_1 = {0.623489801859, -0.781831482468};
                static constexpr T r_2 = {-0.222520933956, -0.974927912182};
                static constexpr T r_3 = {-0.900968867902, -0.433883739118};
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                const T x_3 = forward ? in[3] : conj(in[3]);
                const T x_4 = forward ? in[4] : conj(in[4]);
                const T x_5 = forward ? in[5] : conj(in[5]);
                const T x_6 = forward ? in[6] : conj(in[6]);
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

        template<bool forward>
        struct DFT<8, forward> {
            static void Init() {}
            static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
                static constexpr float c = 0.707106781187f;
                static constexpr T eighth_root = {c,c};
                static constexpr auto forth_root = [](T x){
                    return T{x[1], -x[0]};
                };
                const T x_0 = forward ? in[0] : conj(in[0]);
                const T x_1 = forward ? in[1] : conj(in[1]);
                const T x_2 = forward ? in[2] : conj(in[2]);
                const T x_3 = forward ? in[3] : conj(in[3]);
                const T x_4 = forward ? in[4] : conj(in[4]);
                const T x_5 = forward ? in[5] : conj(in[5]);
                const T x_6 = forward ? in[6] : conj(in[6]);
                const T x_7 = forward ? in[7] : conj(in[7]);
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

        template<unsigned int N, unsigned int S, bool forward>
        static void transpose(T *__restrict__ in, T *__restrict__ out) noexcept {
            constexpr unsigned int A = FFTFactorGen<T>::get_best_factor(N);
            constexpr unsigned int B = N / A;

            // Transpose result back in order
            for (unsigned int i = 0; i < A; i++) {
                if constexpr (B == 1) {
                    // Normalize
                    // FFTW does not normalize, so this division should be commented
                    //  out for appropriate performance comparisons.
                
                    // If the fft is normalized, the ifft does not need to be normalized
                    //  to keep the ifft the exact inverse operation of the fft.
                    out[S * i] = forward ? (in[i] / static_cast<float>(A * S)) : conj(in[i]);
                } else {
                    transpose<B, A * S, forward>(in + (B * i), out + (S * i));
                }
            }
            return;
        }

        template<unsigned int L, unsigned int N>
        struct FFTRecurseLayer {
            static constexpr unsigned int A = FFTFactorGen<T>::get_best_factor(N);
            static constexpr unsigned int B = N / A;
            static constexpr bool base_case = (B == 1);

            static void Init() {
                if constexpr (!base_case) {
                    // Ensure the coeffs are calculated
                    volatile T* coefs = FFTPlan<T>::get_twiddle_factors_by_angle<A, B>();
                    DFT<A, true>::Init();

                    // As well as the next layer's coeffs
                    FFTRecurseLayer<L, B>::Init();
                }
            }

            static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
                // Do the A sized FFT on each bin
                for (unsigned int i = 0; i < L/A; i++) {
                    DFT<A, true>::execute(in + (i * A), out + (i * A));
                }

                if constexpr (!base_case) {
                    // Cooley Tukey twiddle factors
                    T *twiddle_factors = std::assume_aligned<MY_CPU_LOAD_SIZE>(
                        FFTPlan<T>::get_twiddle_factors_by_angle<A, B>());
                    
                    /* This nested loop (1) transposes the data around the second radix B
                     *  to create A bins of size B (2) multiplies the transposed data
                     *  by this layer's twiddle factors (3) transposes the result around 
                     *  what is technically the next layer's A (used as C here) to create 
                     *  D bins (where D is the next layer's B) of size C.
                     *
                     * .All 3 operations are fused to (1) minimize the amount of intermediate
                     *  workspace memory required for the operation (2) minimize overall 
                     *  data movment (CPU IO) (3) structure the operations in a way that is
                     *  easier to tile as a comprehensive unit, if required.
                     *  
                     *  Currently the loops are structured so all reads are consecutive from
                     *  their buffers, but the operations are not tiled.
                     */
                    constexpr unsigned int C = FFTRecurseLayer<L, B>::A;
                    constexpr unsigned int D = FFTRecurseLayer<L, B>::B;
                    for (unsigned int i = 0; i < L/N; i++) {
                        T * const temp_in  = in  + i * N;
                        T * const temp_out = out + i * N;
                        for (unsigned int l = 0; l < C; l++) {
                            for (unsigned int k = 0; k < D; k++) {
                                for (unsigned int j = 0; j < A; j++) {
                                    const unsigned int read_index  = j     + k * A + l * D * A;
                                    const unsigned int write_index = j * B + k * C + l;
                                    temp_in[write_index] =
                                        m(temp_out[read_index], twiddle_factors[read_index]);
                                }
                            }
                        }
                    } 

                    /* Recursively do the B sized FFTs on each bin, which will execute using
                     *  this recursive layer with the next layer's N as this layer's B.
                     *  
                     * Note this layer's B (the next layer's N) is always a factor of this
                     *  layer's N, so inductively N will be decreasing until B is 1. When B
                     *  is 1 we hit our base case and return. This guarantees not only that
                     *  the recursive calls will terminate at runtime, but also that the
                     *  recursive template generation will terminate at compile time.
                     */
                    FFTRecurseLayer<L, B>::fft_recurse(in, out);
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
