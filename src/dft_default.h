#pragma once

  template <unsigned int L, unsigned int F, bool forward, SIMD_TYPE simd>
  struct DFTLayer {
    static void Init() {
      volatile T *coefs = FFTPlan<T, simd>::template get_dft_matrix_by_angle<F>();
      (void)coefs;
    }

    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / F; i++) {
        dft(in + scalar_count(i * F), out + scalar_count(i * F));
      }
    }

    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static T *dft_matrix = static_cast<T *>(__builtin_assume_aligned(
          FFTPlan<T, simd>::template get_dft_matrix_by_angle<F>(),
          MY_CPU_LOAD_SIZE));

      for (unsigned int i = 0; i < F; i++) {
        T sum_real = 0;
        T sum_imag = 0;
        for (unsigned int j = 0; j < F; j++) {
          T in_real, in_imag;
          T mat_real, mat_imag;
          T prod_real, prod_imag;
          load_fft_input<forward>(in, j, in_real, in_imag);
          load(dft_matrix, i * F + j, mat_real, mat_imag);
          multiply(in_real, in_imag, mat_real, mat_imag, prod_real, prod_imag);
          sum_real += prod_real;
          sum_imag += prod_imag;
        }
        store(out, i, sum_real, sum_imag);
      }
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 1, forward, simd> {
    static void Init() {
      static_assert(L != L, "DFTLayer<1> should never be instantiated.");
    }
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      static_assert(L != L, "DFTLayer<1> should never be instantiated.");
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 2, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 2; i++) {
        dft(in + scalar_count(i * 2), out + scalar_count(i * 2));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      T x_0_real, x_0_imag, x_1_real, x_1_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      store(out, 0, x_0_real + x_1_real, x_0_imag + x_1_imag);
      store(out, 1, x_0_real - x_1_real, x_0_imag - x_1_imag);
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 3, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 3; i++) {
        dft(in + scalar_count(i * 3), out + scalar_count(i * 3));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      constexpr T third_root_real = static_cast<T>(-0.5);
      constexpr T third_root_imag = static_cast<T>(-0.866025403784);
      T x_0_real, x_0_imag, x_1_real, x_1_imag, x_2_real, x_2_imag;
      T t1_real, t1_imag, t2_real, t2_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      load_fft_input<forward>(in, 2, x_2_real, x_2_imag);
      store(out, 0, x_0_real + x_1_real + x_2_real,
            x_0_imag + x_1_imag + x_2_imag);
      multiply(x_1_real, x_1_imag, third_root_real, third_root_imag, t1_real,
               t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, third_root_real, third_root_imag,
                        t2_real, t2_imag);
      store(out, 1, x_0_real + t1_real + t2_real,
            x_0_imag + t1_imag + t2_imag);
      multiply_conj_rhs(x_1_real, x_1_imag, third_root_real, third_root_imag,
                        t1_real, t1_imag);
      multiply(x_2_real, x_2_imag, third_root_real, third_root_imag, t2_real,
               t2_imag);
      store(out, 2, x_0_real + t1_real + t2_real,
            x_0_imag + t1_imag + t2_imag);
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 4, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 4; i++) {
        dft(in + scalar_count(i * 4), out + scalar_count(i * 4));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      T x_0_real, x_0_imag, x_1_real, x_1_imag, x_2_real, x_2_imag, x_3_real,
          x_3_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      load_fft_input<forward>(in, 2, x_2_real, x_2_imag);
      load_fft_input<forward>(in, 3, x_3_real, x_3_imag);
      {
        const T a_real = x_0_real + x_2_real;
        const T a_imag = x_0_imag + x_2_imag;
        const T b_real = x_1_real + x_3_real;
        const T b_imag = x_1_imag + x_3_imag;
        store(out, 0, a_real + b_real, a_imag + b_imag);
        store(out, 2, a_real - b_real, a_imag - b_imag);
      }
      {
        const T a_real = x_0_real - x_2_real;
        const T a_imag = x_0_imag - x_2_imag;
        const T temp_b_real = x_1_real - x_3_real;
        const T temp_b_imag = x_1_imag - x_3_imag;
        const T b_real = temp_b_imag;
        const T b_imag = -temp_b_real;
        store(out, 1, a_real + b_real, a_imag + b_imag);
        store(out, 3, a_real - b_real, a_imag - b_imag);
      }
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 5, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 5; i++) {
        dft(in + scalar_count(i * 5), out + scalar_count(i * 5));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      constexpr T root_1_real = static_cast<T>(0.309016994375);
      constexpr T root_1_imag = static_cast<T>(-0.951056516295);
      constexpr T root_2_real = static_cast<T>(-0.809016994375);
      constexpr T root_2_imag = static_cast<T>(-0.587785252292);
      T x_0_real, x_0_imag, x_1_real, x_1_imag, x_2_real, x_2_imag, x_3_real,
          x_3_imag, x_4_real, x_4_imag;
      T t1_real, t1_imag, t2_real, t2_imag, t3_real, t3_imag, t4_real,
          t4_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      load_fft_input<forward>(in, 2, x_2_real, x_2_imag);
      load_fft_input<forward>(in, 3, x_3_real, x_3_imag);
      load_fft_input<forward>(in, 4, x_4_real, x_4_imag);
      store(out, 0, x_0_real + x_1_real + x_2_real + x_3_real + x_4_real,
            x_0_imag + x_1_imag + x_2_imag + x_3_imag + x_4_imag);

      multiply(x_1_real, x_1_imag, root_1_real, root_1_imag, t1_real, t1_imag);
      multiply(x_2_real, x_2_imag, root_2_real, root_2_imag, t2_real, t2_imag);
      multiply_conj_rhs(x_3_real, x_3_imag, root_2_real, root_2_imag, t3_real,
                        t3_imag);
      multiply_conj_rhs(x_4_real, x_4_imag, root_1_real, root_1_imag, t4_real,
                        t4_imag);
      store(out, 1, x_0_real + t1_real + t2_real + t3_real + t4_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag);

      multiply(x_1_real, x_1_imag, root_2_real, root_2_imag, t1_real, t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, root_1_real, root_1_imag, t2_real,
                        t2_imag);
      multiply(x_3_real, x_3_imag, root_1_real, root_1_imag, t3_real, t3_imag);
      multiply_conj_rhs(x_4_real, x_4_imag, root_2_real, root_2_imag, t4_real,
                        t4_imag);
      store(out, 2, x_0_real + t1_real + t2_real + t3_real + t4_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag);

      multiply_conj_rhs(x_1_real, x_1_imag, root_2_real, root_2_imag, t1_real,
                        t1_imag);
      multiply(x_2_real, x_2_imag, root_1_real, root_1_imag, t2_real, t2_imag);
      multiply_conj_rhs(x_3_real, x_3_imag, root_1_real, root_1_imag, t3_real,
                        t3_imag);
      multiply(x_4_real, x_4_imag, root_2_real, root_2_imag, t4_real, t4_imag);
      store(out, 3, x_0_real + t1_real + t2_real + t3_real + t4_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag);

      multiply_conj_rhs(x_1_real, x_1_imag, root_1_real, root_1_imag, t1_real,
                        t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, root_2_real, root_2_imag, t2_real,
                        t2_imag);
      multiply(x_3_real, x_3_imag, root_2_real, root_2_imag, t3_real, t3_imag);
      multiply(x_4_real, x_4_imag, root_1_real, root_1_imag, t4_real, t4_imag);
      store(out, 4, x_0_real + t1_real + t2_real + t3_real + t4_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag);
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 6, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 6; i++) {
        dft(in + scalar_count(i * 6), out + scalar_count(i * 6));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      constexpr T r_a = static_cast<T>(0.5);
      constexpr T r_b = static_cast<T>(-0.866025403784);
      T x_0_real, x_0_imag, x_1_real, x_1_imag, x_2_real, x_2_imag, x_3_real,
          x_3_imag, x_4_real, x_4_imag, x_5_real, x_5_imag;
      T t1_real, t1_imag, t2_real, t2_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      load_fft_input<forward>(in, 2, x_2_real, x_2_imag);
      load_fft_input<forward>(in, 3, x_3_real, x_3_imag);
      load_fft_input<forward>(in, 4, x_4_real, x_4_imag);
      load_fft_input<forward>(in, 5, x_5_real, x_5_imag);
      {
        const T a_real = x_0_real + x_2_real + x_4_real;
        const T a_imag = x_0_imag + x_2_imag + x_4_imag;
        const T b_real = x_1_real + x_3_real + x_5_real;
        const T b_imag = x_1_imag + x_3_imag + x_5_imag;
        store(out, 0, a_real + b_real, a_imag + b_imag);
        store(out, 3, a_real - b_real, a_imag - b_imag);
      }
      {
        const T a_real = x_0_real - x_3_real;
        const T a_imag = x_0_imag - x_3_imag;
        const T b_real = x_1_real - x_4_real;
        const T b_imag = x_1_imag - x_4_imag;
        const T c_real = x_2_real - x_5_real;
        const T c_imag = x_2_imag - x_5_imag;
        multiply(b_real, b_imag, r_a, r_b, t1_real, t1_imag);
        multiply(c_real, c_imag, -r_a, r_b, t2_real, t2_imag);
        store(out, 1, a_real + t1_real + t2_real, a_imag + t1_imag + t2_imag);
        multiply(b_real, b_imag, r_a, -r_b, t1_real, t1_imag);
        multiply(c_real, c_imag, -r_a, -r_b, t2_real, t2_imag);
        store(out, 5, a_real + t1_real + t2_real, a_imag + t1_imag + t2_imag);
      }
      {
        const T a_real = x_0_real + x_3_real;
        const T a_imag = x_0_imag + x_3_imag;
        const T b_real = x_1_real + x_4_real;
        const T b_imag = x_1_imag + x_4_imag;
        const T c_real = x_2_real + x_5_real;
        const T c_imag = x_2_imag + x_5_imag;
        multiply(b_real, b_imag, -r_a, r_b, t1_real, t1_imag);
        multiply(c_real, c_imag, -r_a, -r_b, t2_real, t2_imag);
        store(out, 2, a_real + t1_real + t2_real, a_imag + t1_imag + t2_imag);
        multiply(b_real, b_imag, -r_a, -r_b, t1_real, t1_imag);
        multiply(c_real, c_imag, -r_a, r_b, t2_real, t2_imag);
        store(out, 4, a_real + t1_real + t2_real, a_imag + t1_imag + t2_imag);
      }
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 7, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 7; i++) {
        dft(in + scalar_count(i * 7), out + scalar_count(i * 7));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      constexpr T r_1_real = static_cast<T>(0.623489801859);
      constexpr T r_1_imag = static_cast<T>(-0.781831482468);
      constexpr T r_2_real = static_cast<T>(-0.222520933956);
      constexpr T r_2_imag = static_cast<T>(-0.974927912182);
      constexpr T r_3_real = static_cast<T>(-0.900968867902);
      constexpr T r_3_imag = static_cast<T>(-0.433883739118);
      T x_0_real, x_0_imag, x_1_real, x_1_imag, x_2_real, x_2_imag, x_3_real,
          x_3_imag, x_4_real, x_4_imag, x_5_real, x_5_imag, x_6_real,
          x_6_imag;
      T t1_real, t1_imag, t2_real, t2_imag, t3_real, t3_imag, t4_real,
          t4_imag, t5_real, t5_imag, t6_real, t6_imag;
      load_fft_input<forward>(in, 0, x_0_real, x_0_imag);
      load_fft_input<forward>(in, 1, x_1_real, x_1_imag);
      load_fft_input<forward>(in, 2, x_2_real, x_2_imag);
      load_fft_input<forward>(in, 3, x_3_real, x_3_imag);
      load_fft_input<forward>(in, 4, x_4_real, x_4_imag);
      load_fft_input<forward>(in, 5, x_5_real, x_5_imag);
      load_fft_input<forward>(in, 6, x_6_real, x_6_imag);
      store(out, 0,
            x_0_real + x_1_real + x_2_real + x_3_real + x_4_real + x_5_real +
                x_6_real,
            x_0_imag + x_1_imag + x_2_imag + x_3_imag + x_4_imag + x_5_imag +
                x_6_imag);

      multiply(x_1_real, x_1_imag, r_1_real, r_1_imag, t1_real, t1_imag);
      multiply(x_2_real, x_2_imag, r_2_real, r_2_imag, t2_real, t2_imag);
      multiply(x_3_real, x_3_imag, r_3_real, r_3_imag, t3_real, t3_imag);
      multiply_conj_rhs(x_4_real, x_4_imag, r_3_real, r_3_imag, t4_real,
                        t4_imag);
      multiply_conj_rhs(x_5_real, x_5_imag, r_2_real, r_2_imag, t5_real,
                        t5_imag);
      multiply_conj_rhs(x_6_real, x_6_imag, r_1_real, r_1_imag, t6_real,
                        t6_imag);
      store(out, 1,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);

      multiply(x_1_real, x_1_imag, r_2_real, r_2_imag, t1_real, t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, r_3_real, r_3_imag, t2_real,
                        t2_imag);
      multiply_conj_rhs(x_3_real, x_3_imag, r_1_real, r_1_imag, t3_real,
                        t3_imag);
      multiply(x_4_real, x_4_imag, r_1_real, r_1_imag, t4_real, t4_imag);
      multiply(x_5_real, x_5_imag, r_3_real, r_3_imag, t5_real, t5_imag);
      multiply_conj_rhs(x_6_real, x_6_imag, r_2_real, r_2_imag, t6_real,
                        t6_imag);
      store(out, 2,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);

      multiply(x_1_real, x_1_imag, r_3_real, r_3_imag, t1_real, t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, r_1_real, r_1_imag, t2_real,
                        t2_imag);
      multiply(x_3_real, x_3_imag, r_2_real, r_2_imag, t3_real, t3_imag);
      multiply_conj_rhs(x_4_real, x_4_imag, r_2_real, r_2_imag, t4_real,
                        t4_imag);
      multiply(x_5_real, x_5_imag, r_1_real, r_1_imag, t5_real, t5_imag);
      multiply_conj_rhs(x_6_real, x_6_imag, r_3_real, r_3_imag, t6_real,
                        t6_imag);
      store(out, 3,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);

      multiply_conj_rhs(x_1_real, x_1_imag, r_3_real, r_3_imag, t1_real,
                        t1_imag);
      multiply(x_2_real, x_2_imag, r_1_real, r_1_imag, t2_real, t2_imag);
      multiply_conj_rhs(x_3_real, x_3_imag, r_2_real, r_2_imag, t3_real,
                        t3_imag);
      multiply(x_4_real, x_4_imag, r_2_real, r_2_imag, t4_real, t4_imag);
      multiply_conj_rhs(x_5_real, x_5_imag, r_1_real, r_1_imag, t5_real,
                        t5_imag);
      multiply(x_6_real, x_6_imag, r_3_real, r_3_imag, t6_real, t6_imag);
      store(out, 4,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);

      multiply_conj_rhs(x_1_real, x_1_imag, r_2_real, r_2_imag, t1_real,
                        t1_imag);
      multiply(x_2_real, x_2_imag, r_3_real, r_3_imag, t2_real, t2_imag);
      multiply(x_3_real, x_3_imag, r_1_real, r_1_imag, t3_real, t3_imag);
      multiply_conj_rhs(x_4_real, x_4_imag, r_1_real, r_1_imag, t4_real,
                        t4_imag);
      multiply_conj_rhs(x_5_real, x_5_imag, r_3_real, r_3_imag, t5_real,
                        t5_imag);
      multiply(x_6_real, x_6_imag, r_2_real, r_2_imag, t6_real, t6_imag);
      store(out, 5,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);

      multiply_conj_rhs(x_1_real, x_1_imag, r_1_real, r_1_imag, t1_real,
                        t1_imag);
      multiply_conj_rhs(x_2_real, x_2_imag, r_2_real, r_2_imag, t2_real,
                        t2_imag);
      multiply_conj_rhs(x_3_real, x_3_imag, r_3_real, r_3_imag, t3_real,
                        t3_imag);
      multiply(x_4_real, x_4_imag, r_3_real, r_3_imag, t4_real, t4_imag);
      multiply(x_5_real, x_5_imag, r_2_real, r_2_imag, t5_real, t5_imag);
      multiply(x_6_real, x_6_imag, r_1_real, r_1_imag, t6_real, t6_imag);
      store(out, 6,
            x_0_real + t1_real + t2_real + t3_real + t4_real + t5_real +
                t6_real,
            x_0_imag + t1_imag + t2_imag + t3_imag + t4_imag + t5_imag +
                t6_imag);
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 8, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 8; i++) {
        dft(in + scalar_count(i * 8), out + scalar_count(i * 8));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      constexpr T c = static_cast<T>(0.707106781187);
      T a_0_0_real, a_0_0_imag, a_0_1_real, a_0_1_imag, a_0_2_real,
          a_0_2_imag, a_0_3_real, a_0_3_imag;
      T a_1_0_real, a_1_0_imag, a_1_1_real, a_1_1_imag, a_1_2_real,
          a_1_2_imag, a_1_3_real, a_1_3_imag;
      load_fft_input<forward>(in, 0, a_0_0_real, a_0_0_imag);
      load_fft_input<forward>(in, 1, a_0_1_real, a_0_1_imag);
      load_fft_input<forward>(in, 2, a_0_2_real, a_0_2_imag);
      load_fft_input<forward>(in, 3, a_0_3_real, a_0_3_imag);
      load_fft_input<forward>(in, 4, a_1_0_real, a_1_0_imag);
      load_fft_input<forward>(in, 5, a_1_1_real, a_1_1_imag);
      load_fft_input<forward>(in, 6, a_1_2_real, a_1_2_imag);
      load_fft_input<forward>(in, 7, a_1_3_real, a_1_3_imag);

      const T b_0_0_real = a_0_0_real + a_1_0_real;
      const T b_0_0_imag = a_0_0_imag + a_1_0_imag;
      const T b_0_1_real = a_0_1_real + a_1_1_real;
      const T b_0_1_imag = a_0_1_imag + a_1_1_imag;
      const T b_0_2_real = a_0_2_real + a_1_2_real;
      const T b_0_2_imag = a_0_2_imag + a_1_2_imag;
      const T b_0_3_real = a_0_3_real + a_1_3_real;
      const T b_0_3_imag = a_0_3_imag + a_1_3_imag;

      const T b_1_0_real = a_0_0_real - a_1_0_real;
      const T b_1_0_imag = a_0_0_imag - a_1_0_imag;
      const T b_1_1_real = a_0_1_real - a_1_1_real;
      const T b_1_1_imag = a_0_1_imag - a_1_1_imag;
      const T b_1_2_real = a_0_2_real - a_1_2_real;
      const T b_1_2_imag = a_0_2_imag - a_1_2_imag;
      const T b_1_3_real = a_0_3_real - a_1_3_real;
      const T b_1_3_imag = a_0_3_imag - a_1_3_imag;

      const T c_0_0_real = b_0_0_real + b_0_2_real;
      const T c_0_0_imag = b_0_0_imag + b_0_2_imag;
      const T c_0_1_real = b_0_0_real - b_0_2_real;
      const T c_0_1_imag = b_0_0_imag - b_0_2_imag;
      const T c_0_2_real = b_0_1_real + b_0_3_real;
      const T c_0_2_imag = b_0_1_imag + b_0_3_imag;
      const T c_0_3_real = b_0_1_real - b_0_3_real;
      const T c_0_3_imag = b_0_1_imag - b_0_3_imag;

      T rot_real, rot_imag;
      rotate_by_pos_i(b_1_2_real, b_1_2_imag, rot_real, rot_imag);
      const T c_1_0_real = b_1_0_real + rot_real;
      const T c_1_0_imag = b_1_0_imag + rot_imag;
      const T c_1_1_real = b_1_0_real - rot_real;
      const T c_1_1_imag = b_1_0_imag - rot_imag;
      rotate_by_pos_i(b_1_3_real, b_1_3_imag, rot_real, rot_imag);
      const T c_1_2_real = b_1_1_real + rot_real;
      const T c_1_2_imag = b_1_1_imag + rot_imag;
      const T c_1_3_real = b_1_1_real - rot_real;
      const T c_1_3_imag = b_1_1_imag - rot_imag;

      T t_real, t_imag;
      store(out, 0, c_0_0_real + c_0_2_real, c_0_0_imag + c_0_2_imag);
      multiply(c_1_2_real, c_1_2_imag, c, -c, t_real, t_imag);
      store(out, 1, c_1_0_real + t_real, c_1_0_imag + t_imag);
      rotate_by_pos_i(c_0_3_real, c_0_3_imag, t_real, t_imag);
      store(out, 2, c_0_1_real + t_real, c_0_1_imag + t_imag);
      multiply(c_1_3_real, c_1_3_imag, -c, -c, t_real, t_imag);
      store(out, 3, c_1_1_real + t_real, c_1_1_imag + t_imag);

      store(out, 4, c_0_0_real - c_0_2_real, c_0_0_imag - c_0_2_imag);
      multiply(c_1_2_real, c_1_2_imag, c, -c, t_real, t_imag);
      store(out, 5, c_1_0_real - t_real, c_1_0_imag - t_imag);
      rotate_by_pos_i(c_0_3_real, c_0_3_imag, t_real, t_imag);
      store(out, 6, c_0_1_real - t_real, c_0_1_imag - t_imag);
      multiply(c_1_3_real, c_1_3_imag, -c, -c, t_real, t_imag);
      store(out, 7, c_1_1_real - t_real, c_1_1_imag - t_imag);
    }
  };
