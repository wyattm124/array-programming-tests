#pragma once

#include <cmath>
#include <cstdlib>

#include "prime_factor.hpp"

// TODO: Conditional include on ARCH_x86_64
#include <immintrin.h>

/// TODO:
/// (1) - May want base cases for 9 and 10
/// (2) - Implement Rader's algorithm for small prime numbers, especially 11,
///   13, 17 and 19 
/// (3) - Implement Bluestiens's for numbers that do not decompse
///   into any hand written cases 
/// (4) - Implement a top level search to appropriately pad an input to a
///   number if needed and factor it into the best radices

// For LLVM MCA Analysis
#define MCA_START __asm volatile("# LLVM-MCA-BEGIN");
#define MCA_END __asm volatile("# LLVM-MCA-END");

// For Cache aligning
constexpr unsigned int MY_CPU_LOAD_SIZE = 64;
constexpr unsigned int MY_MAX_SIMD_SIZE = 32;
constexpr unsigned int MY_MAX_ALIGNMENT =
    MY_CPU_LOAD_SIZE > MY_MAX_SIMD_SIZE ? MY_CPU_LOAD_SIZE : MY_MAX_SIMD_SIZE;

enum class SIMD_TYPE : unsigned char {
  NONE = 0,
  AVX2 = 1,
  NEON = 2,
};

namespace FFT {

/// XOR sign-flip masks for `DFT_AVX2_8x1` (bit pattern 0x80000000 == IEEE
///   -0.f).
namespace detail_avx2_dft8 {
// Eighth root of unity, real and imaginary part
constexpr float c = 0.707106781187f;

// mask for identity and negation
constexpr float keep = 0.0f;
constexpr float neg = -0.0f;

alignas(32) inline constexpr float conj_all_mask[] = {keep, neg, keep, neg,
                                                      keep, neg, keep, neg};
alignas(32) inline constexpr float c_0_mask[8] = {keep, keep, neg, neg,
                                                  keep, keep, neg, neg};
alignas(32) inline constexpr float c_1_mask[8] = {keep, neg, neg, keep,
                                                  keep, neg, neg, keep};

// Final mix
alignas(32) inline constexpr float orig_mult[8] = {1.0f, 1.0f, c,  c,
                                                   0.0f, 0.0f, -c, -c};
alignas(32) inline constexpr float flip_mult[8] = {0.0f, 0.0f,  c, -c,
                                                   1.0f, -1.0f, c, -c};
} // namespace detail_avx2_dft8

/// Signal generation for generating an FFT input with the corresponding
/// expected output.
constexpr double TwoPI = 2.0 * M_PI;

template <typename T>
constexpr void wave_gen(T *time_domain, T *freq_domain, unsigned int N,
                        unsigned int f = 1, unsigned int phase = 0,
                        unsigned int amp = 1) {
  for (unsigned int i = 0; i < N; i++) {
    const float angle = static_cast<float>(TwoPI) *
                        static_cast<float>((i * f) + phase) /
                        static_cast<float>(N);
    time_domain[i] += T{cosf(angle), sinf(angle)} * static_cast<float>(amp);
  }
  const float angle = static_cast<float>(TwoPI) *
                      (static_cast<float>(phase) / static_cast<float>(N));
  freq_domain[f] += T{cosf(angle), sinf(angle)} * static_cast<float>(amp);
}

template <typename T>
constexpr void wave_gen_lcg(T *time_domain, T *freq_domain, unsigned int N) {
  if (N < 13) {
    if (N > 7)
      wave_gen(time_domain, freq_domain, N, 7, 2, 1);
    if (N > 5)
      wave_gen(time_domain, freq_domain, N, 5, 2, 1);
    if (N > 4)
      wave_gen(time_domain, freq_domain, N, 4, 3, 2);
    if (N > 3)
      wave_gen(time_domain, freq_domain, N, 3, 2, 1);
    wave_gen(time_domain, freq_domain, N, 1, 1, 1);
    return;
  } else {
    for (unsigned int i = 0; i < 13; i++) {
      wave_gen(time_domain, freq_domain, N, ((i + 7) * 3) % N,
               ((i + 5) * 11) % N, ((i + 11) * 13) % N);
    }
  }
}

/// Operations used by FFT
template <typename T> inline float abs(const T &__restrict__ a) noexcept {
  return sqrtf(a[0] * a[0] + a[1] * a[1]);
}

template <typename T> constexpr T conj(const T &__restrict__ a) noexcept {
  return {a[0], -a[1]};
}

template <typename T> constexpr T flipper(const T &__restrict__ a) noexcept {
  return {a[1], a[0]};
}

template <typename T>
inline T mult(const T &__restrict__ a, const T &__restrict__ b) noexcept {
  return {a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
}
template <typename T>
inline T mult_conj(const T &__restrict__ a,
                   const T &__restrict__ b_conj) noexcept {
  return {a[0] * b_conj[0] + a[1] * b_conj[1],
          a[1] * b_conj[0] - a[0] * b_conj[1]};
}

/// FFT factor generation
template <typename T> struct FFTFactorGen {
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

template <typename T, unsigned int N>
void populate_dft_matrix_by_angle(T *mat) {
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      const double angle =
          NegTwoPI * static_cast<double>((i * j) % N) / static_cast<double>(N);
      mat[i * N + j] =
          T{static_cast<float>(cos(angle)), static_cast<float>(sin(angle))};
    }
  }
}

template <typename T, unsigned int A, unsigned int B>
void populate_twiddle_factors_by_angle(T *factors) {
  for (unsigned int j = 0; j < B; j++) {
    for (unsigned int i = 0; i < A; i++) {
      const double angle =
          NegTwoPI * static_cast<double>(i * j) / static_cast<double>(A * B);
      factors[j * A + i] =
          T{static_cast<float>(cos(angle)), static_cast<float>(sin(angle))};
    }
  }
}

template <typename T, SIMD_TYPE plan_simd = SIMD_TYPE::NONE> class FFTPlan {
public:
  // Multiplication op aliases for type parameter
  static constexpr auto m = mult<T>;
  static constexpr auto mc = mult_conj<T>;
  static constexpr SIMD_TYPE dft_simd = plan_simd;

  // Layers need to be initialized recursively to populate
  //  their corresponding coefficient arrays.
  template <unsigned int N> static void Init() {
    FFTRecurseLayer<N, N>::Init();
  }

  template <unsigned int N>
  static void fft(T *__restrict__ in_unaligned,
                  T *__restrict__ out_unaligned) noexcept {
    T *in = static_cast<T *>(
        __builtin_assume_aligned(in_unaligned, MY_MAX_ALIGNMENT));
    T *out = static_cast<T *>(
        __builtin_assume_aligned(out_unaligned, MY_MAX_ALIGNMENT));
    top_level_fft<N, true>(in, out);
    return;
  }

  template <unsigned int N>
  static void ifft(T *__restrict__ in_unaligned,
                   T *__restrict__ out_unaligned) noexcept {
    T *in = static_cast<T *>(
        __builtin_assume_aligned(in_unaligned, MY_MAX_ALIGNMENT));
    T *out = static_cast<T *>(
        __builtin_assume_aligned(out_unaligned, MY_MAX_ALIGNMENT));
    top_level_fft<N, false>(in, out);
    return;
  }

private:
  // switchable forward / backward fft top level function to reduce redundant
  // code
  template <unsigned int N, bool forward>
  static void top_level_fft(T *__restrict__ in, T *__restrict__ out) noexcept {
    // The result is not shifted so the DC component is still at index 0,
    //  and the FFT and IFFT are inverses of each other

    // The Inverse DFT matrix is the same as the DFT matrix but with the
    // corresponding
    //  factors conjugated.

    if constexpr (FFTRecurseLayer<N, N>::base_case) {
      // If we are already at a base case, then we just need to do a
      // straightforward
      //  DFT.
      DFTLayer<N, N, forward, dft_simd>::execute(in, out);

      for (int i = 0; i < N; i++) {
        out[i] = forward ? (out[i] / static_cast<float>(N)) : conj(out[i]);
      }
    } else {
      // TODO : large FFTs overflow the stack, this array may need to be put on
      // the heap.
      alignas(MY_CPU_LOAD_SIZE) T temp[N];
      constexpr unsigned int A = FFTRecurseLayer<N, N>::A;
      constexpr unsigned int B = FFTRecurseLayer<N, N>::B;

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
      FFTRecurseLayer<N, N>::fft_recurse(out, temp);
      transpose<N, 1, forward>(temp, out);
    }

    return;
  }

  template <unsigned int L, unsigned int F, bool forward,
            SIMD_TYPE simd = SIMD_TYPE::NONE>
  struct DFTLayer;

#include "dft_default.h"
#include "dft_AVX2.h"

  template <unsigned int N, unsigned int S, bool forward>
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
        out[S * i] =
            forward ? (in[i] / static_cast<float>(A * S)) : conj(in[i]);
      } else {
        transpose<B, A * S, forward>(in + (B * i), out + (S * i));
      }
    }
    return;
  }

  template <unsigned int L, unsigned int N, unsigned int A, unsigned int B,
            unsigned int C, unsigned int D>
  static void inner_transpose(T *__restrict__ in, T *__restrict__ out,
                              T *__restrict__ twiddle_factors) {
    for (unsigned int i = 0; i < L / N; i++) {
      T *const temp_in = in + i * N;
      T *const temp_out = out + i * N;
      for (unsigned int l = 0; l < C; l++) {
        for (unsigned int k = 0; k < D; k++) {
          for (unsigned int j = 0; j < A; j++) {
            const unsigned int read_index = j + k * A + l * D * A;
            const unsigned int write_index = j * B + k * C + l;
            temp_out[write_index] =
                m(temp_in[read_index], twiddle_factors[read_index]);
          }
        }
      }
    }
  }

  template <unsigned int L, unsigned int N> struct FFTRecurseLayer {
    static constexpr unsigned int A = FFTFactorGen<T>::get_best_factor(N);
    static constexpr unsigned int B = N / A;
    static constexpr bool base_case = (B == 1);

    static void Init() {
      if constexpr (!base_case) {
        // Ensure the coeffs are calculated
        volatile T *coefs =
            FFTPlan<T, plan_simd>::template get_twiddle_factors_by_angle<A, B>();
        DFTLayer<L, A, true, dft_simd>::Init();

        // As well as the next layer's coeffs
        FFTRecurseLayer<L, B>::Init();
      }
    }

    static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
      // Do the A sized FFT on each bin
      DFTLayer<L, A, true, dft_simd>::execute(in, out);

      if constexpr (!base_case) {
        // Cooley Tukey twiddle factors
        T *twiddle_factors = static_cast<T *>(__builtin_assume_aligned(
            FFTPlan<T, plan_simd>::template get_twiddle_factors_by_angle<A, B>(),
            MY_CPU_LOAD_SIZE));

        /* This nested loop (1) transposes the data around the second radix B
         *  to create A bins of size B (2) multiplies the transposed data
         *  by this layer's twiddle factors (3) transposes the result around
         *  what is technically the next layer's A (used as C here) to create
         *  D bins (where D is the next layer's B) of size C.
         *
         * All 3 operations are fused to (1) minimize the amount of intermediate
         *  workspace memory required for the operation (2) minimize overall
         *  data movment (CPU IO) (3) structure the operations in a way that is
         *  easier to tile as a comprehensive unit, if required.
         *
         *  Currently the loops are structured so all reads are consecutive from
         *  their buffers, but the operations are not tiled.
         */

        inner_transpose<L, N, A, B, FFTRecurseLayer<L, B>::A,
                        FFTRecurseLayer<L, B>::B>(out, in, twiddle_factors);

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

  // Coeffs calculated ahead of time trades memory for speed.
  //   Since Coeffs are roots of unity they could also be multiplied
  //   together within loops to generate them as needed, but this
  //   has emperically been found to be not acceptably stable numerically.
  //   It is better to generate them with this angle based method for
  //   numerical stablility.

  // TODO : use a cleaner way to hold coefs at class rather than object level
  template <unsigned int N> static T *get_dft_matrix_by_angle() {
    static T *coefs = [] {
      T *result =
          static_cast<T *>(aligned_alloc(MY_MAX_ALIGNMENT, sizeof(T) * N * N));
      populate_dft_matrix_by_angle<T, N>(result);
      return result;
    }();
    return coefs;
  };

  template <unsigned int A, unsigned int B>
  static T *get_twiddle_factors_by_angle() {
    static T *coefs = [] {
      T *result =
          static_cast<T *>(aligned_alloc(MY_MAX_ALIGNMENT, sizeof(T) * A * B));
      populate_twiddle_factors_by_angle<T, A, B>(result);
      return result;
    }();
    return coefs;
  };
};
} // namespace FFT
