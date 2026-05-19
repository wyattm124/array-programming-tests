#pragma once

  // A full DFT done by simple matrix multiplication with the appropriate DFT
  // matrix.
  //  Specific sizes are specialized and optimized as appropriate.
  template <unsigned int L, unsigned int F, bool forward, SIMD_TYPE simd>
  struct DFTLayer {
    static void Init() {
      // Ensure the coeffs are calculated
      volatile T *coefs =
          FFTPlan<T, simd>::template get_dft_matrix_by_angle<F>();
    }

    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / F; i++) {
        dft(in + i * F, out + i * F);
      }
    }

    // Simply compute the multiplication of the input array (vector) by
    //  by the DFT matrix.
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static T *dft_matrix = static_cast<T *>(__builtin_assume_aligned(
          FFTPlan<T, simd>::template get_dft_matrix_by_angle<F>(),
          MY_CPU_LOAD_SIZE));

      for (unsigned int i = 0; i < F; i++) {
        T temp_ans = {0, 0};
        for (unsigned int j = 0; j < F; j++) {
          temp_ans += m(forward ? in[j] : conj(in[j]), dft_matrix[i * F + j]);
        }
        out[i] = temp_ans;
      }
    }
  };

  // Hand written DFT cases serve as base cases for recursive Cooley Tookey
  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 1, forward, simd> {
    static void Init() {
      // This init should never be called, as this FFTRecurseLayer should never
      // be used
      static_assert(false);
    }
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      // All base cases should avoid the need for this trivial specialization
      static_assert(false);
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 2, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 2; i++) {
        dft(in + i * 2, out + i * 2);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      out[0] = x_0 + x_1;
      out[1] = x_0 - x_1;
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 3, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 3; i++) {
        dft(in + i * 3, out + i * 3);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static constexpr T third_root = {-0.5f, -0.866025403784f};
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      const T &x_2 = forward ? in[2] : conj(in[2]);
      out[0] = x_0 + x_1 + x_2;
      out[1] = x_0 + m(x_1, third_root) + mc(x_2, third_root);
      out[2] = x_0 + mc(x_1, third_root) + m(x_2, third_root);
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 4, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 4; i++) {
        dft(in + 4 * i, out + i * 4);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      const T &x_2 = forward ? in[2] : conj(in[2]);
      const T &x_3 = forward ? in[3] : conj(in[3]);
      { // First do evens
        const T &a = x_0 + x_2;
        const T &b = x_1 + x_3;
        out[0] = a + b;
        out[2] = a - b;
      }
      { // Then do odds
        const T &a = x_0 - x_2;
        const T &temp_b = x_1 - x_3;
        const T &b = {temp_b[1], -temp_b[0]};
        out[1] = a + b;
        out[3] = a - b;
      }
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 5, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 5; i++) {
        dft(in + i * 5, out + i * 5);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static constexpr T root_1 = {0.309016994375, -0.951056516295};
      static constexpr T root_2 = {-0.809016994375, -0.587785252292};
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      const T &x_2 = forward ? in[2] : conj(in[2]);
      const T &x_3 = forward ? in[3] : conj(in[3]);
      const T &x_4 = forward ? in[4] : conj(in[4]);
      out[0] = x_0 + x_1 + x_2 + x_3 + x_4;
      out[1] = x_0 + m(x_1, root_1) + m(x_2, root_2) + mc(x_3, root_2) +
               mc(x_4, root_1);
      out[2] = x_0 + m(x_1, root_2) + mc(x_2, root_1) + m(x_3, root_1) +
               mc(x_4, root_2);
      out[3] = x_0 + mc(x_1, root_2) + m(x_2, root_1) + mc(x_3, root_1) +
               m(x_4, root_2);
      out[4] = x_0 + mc(x_1, root_1) + mc(x_2, root_2) + m(x_3, root_2) +
               m(x_4, root_1);
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 6, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 6; i++) {
        dft(in + 6 * i, out + 6 * i);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static constexpr float r_a = 0.5f;
      static constexpr float r_b = -0.866025403784f;
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      const T &x_2 = forward ? in[2] : conj(in[2]);
      const T &x_3 = forward ? in[3] : conj(in[3]);
      const T &x_4 = forward ? in[4] : conj(in[4]);
      const T &x_5 = forward ? in[5] : conj(in[5]);
      {
        const T &a = x_0 + x_2 + x_4;
        const T &b = x_1 + x_3 + x_5;
        out[0] = a + b;
        out[3] = a - b;
      }
      {
        const T &a = x_0 - x_3;
        const T &b = x_1 - x_4;
        const T &c = x_2 - x_5;
        out[1] = a + m(b, T{r_a, r_b}) + m(c, T{-r_a, r_b});
        out[5] = a + m(b, T{r_a, -r_b}) + m(c, T{-r_a, -r_b});
      }
      {
        const T &a = x_0 + x_3;
        const T &b = x_1 + x_4;
        const T &c = x_2 + x_5;
        out[2] = a + m(b, T{-r_a, r_b}) + m(c, T{-r_a, -r_b});
        out[4] = a + m(b, T{-r_a, -r_b}) + m(c, T{-r_a, r_b});
      }
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 7, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 7; i++) {
        dft(in + 7 * i, out + 7 * i);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static constexpr T r_1 = {0.623489801859, -0.781831482468};
      static constexpr T r_2 = {-0.222520933956, -0.974927912182};
      static constexpr T r_3 = {-0.900968867902, -0.433883739118};
      const T &x_0 = forward ? in[0] : conj(in[0]);
      const T &x_1 = forward ? in[1] : conj(in[1]);
      const T &x_2 = forward ? in[2] : conj(in[2]);
      const T &x_3 = forward ? in[3] : conj(in[3]);
      const T &x_4 = forward ? in[4] : conj(in[4]);
      const T &x_5 = forward ? in[5] : conj(in[5]);
      const T &x_6 = forward ? in[6] : conj(in[6]);
      out[0] = x_0 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6;
      out[1] = x_0 + m(x_1, r_1) + m(x_2, r_2) + m(x_3, r_3) + mc(x_4, r_3) +
               mc(x_5, r_2) + mc(x_6, r_1);
      out[2] = x_0 + m(x_1, r_2) + mc(x_2, r_3) + mc(x_3, r_1) + m(x_4, r_1) +
               m(x_5, r_3) + mc(x_6, r_2);
      out[3] = x_0 + m(x_1, r_3) + mc(x_2, r_1) + m(x_3, r_2) + mc(x_4, r_2) +
               m(x_5, r_1) + mc(x_6, r_3);
      out[4] = x_0 + mc(x_1, r_3) + m(x_2, r_1) + mc(x_3, r_2) + m(x_4, r_2) +
               mc(x_5, r_1) + m(x_6, r_3);
      out[5] = x_0 + mc(x_1, r_2) + m(x_2, r_3) + m(x_3, r_1) + mc(x_4, r_1) +
               mc(x_5, r_3) + m(x_6, r_2);
      out[6] = x_0 + mc(x_1, r_1) + mc(x_2, r_2) + mc(x_3, r_3) + m(x_4, r_3) +
               m(x_5, r_2) + m(x_6, r_1);
      return;
    }
  };

  template <unsigned int L, bool forward, SIMD_TYPE simd>
  struct DFTLayer<L, 8, forward, simd> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 8; i++) {
        dft(in + 8 * i, out + 8 * i);
      }
      return;
    }
    static void dft(T *__restrict__ in_unaligned,
                    T *__restrict__ out_unaligned) noexcept {
      T *in = static_cast<T *>(
          __builtin_assume_aligned(in_unaligned, MY_MAX_ALIGNMENT));
      T *out = static_cast<T *>(
          __builtin_assume_aligned(out_unaligned, MY_MAX_ALIGNMENT));
      static constexpr float c = 0.707106781187f;
      static constexpr T eighth_root = {c, c};
      static constexpr T eighth_root_conj = {c, -c};
      static constexpr T neg_eighth_root = {-c, -c};
      static constexpr auto fourth_root = [](T x) { return T{x[1], -x[0]}; };

      const T &a_0_0 = forward ? in[0] : conj(in[0]);
      const T &a_0_1 = forward ? in[1] : conj(in[1]);
      const T &a_0_2 = forward ? in[2] : conj(in[2]);
      const T &a_0_3 = forward ? in[3] : conj(in[3]);

      const T &a_1_0 = forward ? in[4] : conj(in[4]);
      const T &a_1_1 = forward ? in[5] : conj(in[5]);
      const T &a_1_2 = forward ? in[6] : conj(in[6]);
      const T &a_1_3 = forward ? in[7] : conj(in[7]);
      {
        const T &b_0_0 = a_0_0 + a_1_0;
        const T &b_0_1 = a_0_1 + a_1_1;
        const T &b_0_2 = a_0_2 + a_1_2;
        const T &b_0_3 = a_0_3 + a_1_3;

        const T &b_1_0 = a_0_0 - a_1_0;
        const T &b_1_1 = a_0_1 - a_1_1;
        const T &b_1_2 = a_0_2 - a_1_2;
        const T &b_1_3 = a_0_3 - a_1_3;
        {
          const T &c_0_0 = b_0_0 + b_0_2;
          const T &c_0_1 = b_0_0 - b_0_2;
          const T &c_0_2 = b_0_1 + b_0_3;
          const T &c_0_3 = b_0_1 - b_0_3;

          const T &c_1_0 = b_1_0 + fourth_root(b_1_2);
          const T &c_1_1 = b_1_0 - fourth_root(b_1_2);
          const T &c_1_2 = b_1_1 + fourth_root(b_1_3);
          const T &c_1_3 = b_1_1 - fourth_root(b_1_3);
          {
            out[0] = c_0_0 + c_0_2;
            out[1] = c_1_0 + m(c_1_2, eighth_root_conj);
            out[2] = c_0_1 + fourth_root(c_0_3);
            out[3] = c_1_1 + m(c_1_3, neg_eighth_root);

            out[4] = c_0_0 - c_0_2;
            out[5] = c_1_0 - m(c_1_2, eighth_root_conj);
            out[6] = c_0_1 - fourth_root(c_0_3);
            out[7] = c_1_1 - m(c_1_3, neg_eighth_root);
          }
        }
      }
      return;
    }
  };
