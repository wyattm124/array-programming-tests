#pragma once 

  template <unsigned int L>
  struct DFTLayer<L, 2, true, SIMD_TYPE::AVX2> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 2; i++) {
        dft(in + i * 2, out + i * 2);
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      static const __m256 neg_mask =
          _mm256_set_ps(0.0f, 0.0f, -0.0f, -0.0f, 0.0f, 0.0f, -0.0f, -0.0f);
      static const __m256i flipper_indices =
          _mm256_set_epi32(2, 3, 0, 1, 6, 7, 4, 5);

      const __m256 x = _mm256_loadu_ps(reinterpret_cast<float *>(in));
      const __m256 res =
          _mm256_add_ps(_mm256_xor_ps(x, neg_mask),
                        _mm256_permutevar8x32_ps(x, flipper_indices));
      _mm256_storeu_ps(reinterpret_cast<float *>(out), res);
      return;
    }
  }; 

  template <unsigned int L, bool forward>
  struct DFTLayer<L, 8, forward, SIMD_TYPE::AVX2> {
    static void Init() {}
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 8; i++) {
        dft(in + 8 * i, out + 8 * i);
      }
      return;
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      MCA_START;
      __m256 a_0 = _mm256_load_ps(reinterpret_cast<float *>(in));
      __m256 a_1 = _mm256_load_ps(reinterpret_cast<float *>(in + 4));

      if constexpr (!forward) {
        a_0 =
            _mm256_xor_ps(a_0, _mm256_load_ps(detail_avx2_dft8::conj_all_mask));
        a_1 =
            _mm256_xor_ps(a_1, _mm256_load_ps(detail_avx2_dft8::conj_all_mask));
      }

      {
        const __m256 b_0 = _mm256_add_ps(a_0, a_1);
        const __m256 b_1 = _mm256_sub_ps(a_0, a_1);

        {
          const __m256 c_0 = _mm256_add_ps(
              _mm256_permute4x64_pd(_mm256_castps_pd(b_0), 0x50),
              _mm256_xor_ps(_mm256_permute4x64_pd(_mm256_castps_pd(b_0), 0xFA),
                            _mm256_load_ps(detail_avx2_dft8::c_0_mask)));

          const __m256 c_1 = _mm256_add_ps(
              _mm256_permute4x64_pd(_mm256_castps_pd(b_1), 0x50),
              _mm256_xor_ps(
                  _mm256_permute_ps(
                      _mm256_permute4x64_pd(_mm256_castps_pd(b_1), 0xFA), 0xB1),
                  _mm256_load_ps(detail_avx2_dft8::c_1_mask)));

          {
            const __m256 temp_d_0 = _mm256_unpacklo_pd(c_0, c_1);
            const __m256 temp_temp_d_1 = _mm256_unpackhi_pd(c_0, c_1);
            const __m256 d_0 =
                _mm256_permute2f128_pd(temp_d_0, temp_temp_d_1, 0x20);
            const __m256 temp_d_1 =
                _mm256_permute2f128_pd(temp_d_0, temp_temp_d_1, 0x31);

            const __m256 d_1 = _mm256_fmadd_ps(
                temp_d_1, _mm256_load_ps(detail_avx2_dft8::orig_mult),
                _mm256_mul_ps(_mm256_permute_ps(temp_d_1, 0xB1),
                              _mm256_load_ps(detail_avx2_dft8::flip_mult)));

            _mm256_store_ps(reinterpret_cast<float *>(out),
                            _mm256_add_ps(d_0, d_1));
            _mm256_store_ps(reinterpret_cast<float *>(out + 4),
                            _mm256_sub_ps(d_0, d_1));
          }
        }
      }
      MCA_END;
      return;
    }
  };
