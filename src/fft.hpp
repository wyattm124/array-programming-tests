#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

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

constexpr unsigned int COMPLEX_STRIDE = 2;

enum class SIMD_TYPE : unsigned char {
  NONE = 0,
  AVX2 = 1,
  NEON = 2,
};

namespace FFT {

/// XOR sign-flip masks for `DFT_AVX2_8x1` (bit pattern 0x80000000 == IEEE
///   -0.f).
namespace detail_avx2_dft8 {
constexpr float c = 0.707106781187f;
constexpr float keep = 0.0f;
constexpr float neg = -0.0f;

alignas(32) inline constexpr float conj_all_mask[] = {keep, neg, keep, neg,
                                                      keep, neg, keep, neg};
alignas(32) inline constexpr float c_0_mask[8] = {keep, keep, neg, neg,
                                                  keep, keep, neg, neg};
alignas(32) inline constexpr float c_1_mask[8] = {keep, neg, neg, keep,
                                                  keep, neg, neg, keep};

alignas(32) inline constexpr float orig_mult[8] = {1.0f, 1.0f, c,  c,
                                                   0.0f, 0.0f, -c, -c};
alignas(32) inline constexpr float flip_mult[8] = {0.0f, 0.0f,  c, -c,
                                                   1.0f, -1.0f, c, -c};
} // namespace detail_avx2_dft8

constexpr double TwoPI = 2.0 * M_PI;
constexpr double NegTwoPI = -TwoPI;

template <typename T>
constexpr unsigned int real_index(unsigned int complex_index) noexcept {
  return complex_index * COMPLEX_STRIDE;
}

template <typename T>
constexpr unsigned int imag_index(unsigned int complex_index) noexcept {
  return real_index<T>(complex_index) + 1;
}

template <typename T>
constexpr T load_real_interleaved(const T *data, unsigned int index) noexcept {
  return data[real_index<T>(index)];
}

template <typename T>
constexpr T load_imag_interleaved(const T *data, unsigned int index) noexcept {
  return data[imag_index<T>(index)];
}

template <typename T>
constexpr void store_interleaved(T *data, unsigned int index, T real,
                                 T imag) noexcept {
  data[real_index<T>(index)] = real;
  data[imag_index<T>(index)] = imag;
}

template <typename T>
inline auto abs(const T &__restrict__ a) noexcept
    -> std::remove_cvref_t<decltype(a[0])> {
  using Scalar = std::remove_cvref_t<decltype(a[0])>;
  using std::sqrt;
  return static_cast<Scalar>(sqrt(a[0] * a[0] + a[1] * a[1]));
}

template <typename T>
  requires std::is_floating_point_v<T>
constexpr void wave_gen(T *time_domain, T *freq_domain, unsigned int N,
                        unsigned int f = 1, unsigned int phase = 0,
                        unsigned int amp = 1) {
  for (unsigned int i = 0; i < N; i++) {
    const T angle = static_cast<T>(TwoPI) *
                    static_cast<T>((i * f) + phase) / static_cast<T>(N);
    const T wave_real = static_cast<T>(std::cos(angle)) * static_cast<T>(amp);
    const T wave_imag = static_cast<T>(std::sin(angle)) * static_cast<T>(amp);
    store_interleaved(
        time_domain, i,
        load_real_interleaved(time_domain, i) + wave_real,
        load_imag_interleaved(time_domain, i) + wave_imag);
  }

  const T angle = static_cast<T>(TwoPI) *
                  (static_cast<T>(phase) / static_cast<T>(N));
  const T wave_real = static_cast<T>(std::cos(angle)) * static_cast<T>(amp);
  const T wave_imag = static_cast<T>(std::sin(angle)) * static_cast<T>(amp);
  store_interleaved(freq_domain, f,
                    load_real_interleaved(freq_domain, f) + wave_real,
                    load_imag_interleaved(freq_domain, f) + wave_imag);
}

template <typename T>
  requires (!std::is_floating_point_v<T>)
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
  }

  for (unsigned int i = 0; i < 13; i++) {
    wave_gen(time_domain, freq_domain, N, ((i + 7) * 3) % N,
             ((i + 5) * 11) % N, ((i + 11) * 13) % N);
  }
}

template <typename T> struct FFTFactorGen {
  static constexpr unsigned int factors[] = {8, 4, 6};

  static constexpr unsigned int get_best_factor(unsigned int N) {
    for (auto factor : factors) {
      if (N % factor == 0)
        return factor;
    }
    return prime_factor::get_prime_factor(N);
  }
};

template <typename T, unsigned int N>
void populate_dft_matrix_by_angle(T *mat) {
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      const double angle =
          NegTwoPI * static_cast<double>((i * j) % N) / static_cast<double>(N);
      store_interleaved(mat, i * N + j, static_cast<T>(std::cos(angle)),
                        static_cast<T>(std::sin(angle)));
    }
  }
}

template <typename T, unsigned int A, unsigned int B>
void populate_twiddle_factors_by_angle(T *factors) {
  for (unsigned int j = 0; j < B; j++) {
    for (unsigned int i = 0; i < A; i++) {
      const double angle =
          NegTwoPI * static_cast<double>(i * j) / static_cast<double>(A * B);
      store_interleaved(factors, j * A + i, static_cast<T>(std::cos(angle)),
                        static_cast<T>(std::sin(angle)));
    }
  }
}

template <typename T, SIMD_TYPE plan_simd = SIMD_TYPE::NONE> class FFTPlan {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                "FFTPlan<T> expects T to be float or double.");

public:
  static constexpr SIMD_TYPE dft_simd = plan_simd;

  template <unsigned int N> static void Init() { FFTRecurseLayer<N, N>::Init(); }

  template <unsigned int N>
  static void fft(T *__restrict__ in_unaligned,
                  T *__restrict__ out_unaligned) noexcept {
    T *in = static_cast<T *>(
        __builtin_assume_aligned(in_unaligned, MY_MAX_ALIGNMENT));
    T *out = static_cast<T *>(
        __builtin_assume_aligned(out_unaligned, MY_MAX_ALIGNMENT));
    top_level_fft<N, true>(in, out);
  }

  template <unsigned int N>
  static void ifft(T *__restrict__ in_unaligned,
                   T *__restrict__ out_unaligned) noexcept {
    T *in = static_cast<T *>(
        __builtin_assume_aligned(in_unaligned, MY_MAX_ALIGNMENT));
    T *out = static_cast<T *>(
        __builtin_assume_aligned(out_unaligned, MY_MAX_ALIGNMENT));
    top_level_fft<N, false>(in, out);
  }

private:
  static constexpr unsigned int scalar_count(unsigned int complex_count) {
    return complex_count * COMPLEX_STRIDE;
  }

  static T real(const T *__restrict__ data, unsigned int index) noexcept {
    return data[real_index<T>(index)];
  }

  static T imag(const T *__restrict__ data, unsigned int index) noexcept {
    return data[imag_index<T>(index)];
  }

  static void load(const T *__restrict__ data, unsigned int index, T &real_part,
                   T &imag_part) noexcept {
    real_part = real(data, index);
    imag_part = imag(data, index);
  }

  template <bool forward>
  static void load_fft_input(const T *__restrict__ data, unsigned int index,
                             T &real_part, T &imag_part) noexcept {
    load(data, index, real_part, imag_part);
    if constexpr (!forward)
      imag_part = -imag_part;
  }

  static void store(T *__restrict__ data, unsigned int index, T real_part,
                    T imag_part) noexcept {
    data[real_index<T>(index)] = real_part;
    data[imag_index<T>(index)] = imag_part;
  }

  static void multiply(T a_real, T a_imag, T b_real, T b_imag, T &out_real,
                       T &out_imag) noexcept {
    out_real = a_real * b_real - a_imag * b_imag;
    out_imag = a_real * b_imag + a_imag * b_real;
  }

  static void multiply_conj_rhs(T a_real, T a_imag, T b_real, T b_imag,
                                T &out_real, T &out_imag) noexcept {
    out_real = a_real * b_real + a_imag * b_imag;
    out_imag = a_imag * b_real - a_real * b_imag;
  }

  static void rotate_by_pos_i(T in_real, T in_imag, T &out_real,
                              T &out_imag) noexcept {
    out_real = in_imag;
    out_imag = -in_real;
  }

  static T *allocate_aligned(std::size_t scalar_elements) {
    const std::size_t bytes = sizeof(T) * scalar_elements;
    const std::size_t aligned_bytes =
        ((bytes + MY_MAX_ALIGNMENT - 1) / MY_MAX_ALIGNMENT) * MY_MAX_ALIGNMENT;
    return static_cast<T *>(aligned_alloc(MY_MAX_ALIGNMENT, aligned_bytes));
  }

  template <unsigned int N, bool forward>
  static void top_level_fft(T *__restrict__ in, T *__restrict__ out) noexcept {
    if constexpr (FFTRecurseLayer<N, N>::base_case) {
      DFTLayer<N, N, forward, dft_simd>::execute(in, out);

      for (unsigned int i = 0; i < N; i++) {
        T out_real, out_imag;
        load(out, i, out_real, out_imag);
        if constexpr (forward) {
          store(out, i, out_real / static_cast<T>(N),
                out_imag / static_cast<T>(N));
        } else {
          store(out, i, out_real, -out_imag);
        }
      }
    } else {
      alignas(MY_CPU_LOAD_SIZE) T temp[scalar_count(N)];
      constexpr unsigned int A = FFTRecurseLayer<N, N>::A;
      constexpr unsigned int B = FFTRecurseLayer<N, N>::B;

      for (unsigned int i = 0; i < B; i++) {
        for (unsigned int j = 0; j < A; j++) {
          T in_real, in_imag;
          load(in, i + j * B, in_real, in_imag);
          if constexpr (!forward)
            in_imag = -in_imag;
          store(out, i * A + j, in_real, in_imag);
        }
      }
      FFTRecurseLayer<N, N>::fft_recurse(out, temp);
      transpose<N, 1, forward>(temp, out);
    }
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

    for (unsigned int i = 0; i < A; i++) {
      if constexpr (B == 1) {
        T in_real, in_imag;
        load(in, i, in_real, in_imag);
        if constexpr (forward) {
          store(out, S * i, in_real / static_cast<T>(A * S),
                in_imag / static_cast<T>(A * S));
        } else {
          store(out, S * i, in_real, -in_imag);
        }
      } else {
        transpose<B, A * S, forward>(in + scalar_count(B * i),
                                     out + scalar_count(S * i));
      }
    }
  }

  template <unsigned int L, unsigned int N, unsigned int A, unsigned int B,
            unsigned int C, unsigned int D>
  static void inner_transpose(T *__restrict__ in, T *__restrict__ out,
                              T *__restrict__ twiddle_factors) {
    for (unsigned int i = 0; i < L / N; i++) {
      T *const temp_in = in + scalar_count(i * N);
      T *const temp_out = out + scalar_count(i * N);
      for (unsigned int l = 0; l < C; l++) {
        for (unsigned int k = 0; k < D; k++) {
          for (unsigned int j = 0; j < A; j++) {
            const unsigned int read_index = j + k * A + l * D * A;
            const unsigned int write_index = j * B + k * C + l;
            T in_real, in_imag;
            T twiddle_real, twiddle_imag;
            T out_real, out_imag;
            load(temp_in, read_index, in_real, in_imag);
            load(twiddle_factors, read_index, twiddle_real, twiddle_imag);
            multiply(in_real, in_imag, twiddle_real, twiddle_imag, out_real,
                     out_imag);
            store(temp_out, write_index, out_real, out_imag);
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
        volatile T *coefs =
            FFTPlan<T, plan_simd>::template get_twiddle_factors_by_angle<A, B>();
        (void)coefs;
        DFTLayer<L, A, true, dft_simd>::Init();
        FFTRecurseLayer<L, B>::Init();
      }
    }

    static void fft_recurse(T *__restrict__ in, T *__restrict__ out) noexcept {
      DFTLayer<L, A, true, dft_simd>::execute(in, out);

      if constexpr (!base_case) {
        T *twiddle_factors = static_cast<T *>(__builtin_assume_aligned(
            FFTPlan<T, plan_simd>::template get_twiddle_factors_by_angle<A, B>(),
            MY_CPU_LOAD_SIZE));

        inner_transpose<L, N, A, B, FFTRecurseLayer<L, B>::A,
                        FFTRecurseLayer<L, B>::B>(out, in, twiddle_factors);
        FFTRecurseLayer<L, B>::fft_recurse(in, out);
      }
    }
  };

  template <unsigned int N> static T *get_dft_matrix_by_angle() {
    static T *coefs = [] {
      T *result = allocate_aligned(static_cast<std::size_t>(N) * N *
                                   COMPLEX_STRIDE);
      populate_dft_matrix_by_angle<T, N>(result);
      return result;
    }();
    return coefs;
  }

  template <unsigned int A, unsigned int B>
  static T *get_twiddle_factors_by_angle() {
    static T *coefs = [] {
      T *result =
          allocate_aligned(static_cast<std::size_t>(A) * B * COMPLEX_STRIDE);
      populate_twiddle_factors_by_angle<T, A, B>(result);
      return result;
    }();
    return coefs;
  }
};
} // namespace FFT
