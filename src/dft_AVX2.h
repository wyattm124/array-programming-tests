#pragma once

  template <unsigned int L>
  struct DFTLayer<L, 2, true, SIMD_TYPE::AVX2> {
    static void Init() {
      static_assert(std::is_same_v<T, float>,
                    "AVX2 DFTLayer specializations currently support float only.");
    }
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 2; i++) {
        dft(in + scalar_count(i * 2), out + scalar_count(i * 2));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      const __m128 x = _mm_loadu_ps(in);
      const __m128 swapped = _mm_shuffle_ps(x, x, _MM_SHUFFLE(1, 0, 3, 2));
      _mm_storeu_ps(out, _mm_add_ps(_mm_mul_ps(x, _mm_set_ps(-1.f, -1.f, 1.f, 1.f)), swapped));
    }
  };

  template <unsigned int L, bool forward>
  struct DFTLayer<L, 8, forward, SIMD_TYPE::AVX2> {
    static void Init() {
      static_assert(std::is_same_v<T, float>,
                    "AVX2 DFTLayer specializations currently support float only.");
    }
    static void execute(T *__restrict__ in, T *__restrict__ out) noexcept {
      for (unsigned int i = 0; i < L / 8; i++) {
        dft(in + scalar_count(i * 8), out + scalar_count(i * 8));
      }
    }
    static void dft(T *__restrict__ in, T *__restrict__ out) noexcept {
      MCA_START;
      __m256 a_0 = _mm256_load_ps(in);
      __m256 a_1 = _mm256_load_ps(in + 8);

      if constexpr (!forward) {
        a_0 = _mm256_xor_ps(a_0, _mm256_load_ps(detail_avx2_dft8::conj_all_mask));
        a_1 = _mm256_xor_ps(a_1, _mm256_load_ps(detail_avx2_dft8::conj_all_mask));
      }

      const __m256 b_0 = _mm256_add_ps(a_0, a_1);
      const __m256 b_1 = _mm256_sub_ps(a_0, a_1);

      const __m256 c_0 = _mm256_add_ps(
          _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(b_0), 0x50)),
          _mm256_xor_ps(
              _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(b_0), 0xFA)),
              _mm256_load_ps(detail_avx2_dft8::c_0_mask)));

      const __m256 c_1 = _mm256_add_ps(
          _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(b_1), 0x50)),
          _mm256_xor_ps(
              _mm256_permute_ps(
                  _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(b_1), 0xFA)),
                  0xB1),
              _mm256_load_ps(detail_avx2_dft8::c_1_mask)));

      const __m256d temp_d_0 = _mm256_unpacklo_pd(_mm256_castps_pd(c_0), _mm256_castps_pd(c_1));
      const __m256d temp_temp_d_1 = _mm256_unpackhi_pd(_mm256_castps_pd(c_0), _mm256_castps_pd(c_1));
      const __m256 d_0 = _mm256_castpd_ps(_mm256_permute2f128_pd(temp_d_0, temp_temp_d_1, 0x20));
      const __m256 temp_d_1 = _mm256_castpd_ps(_mm256_permute2f128_pd(temp_d_0, temp_temp_d_1, 0x31));

      const __m256 d_1 = _mm256_fmadd_ps(
          temp_d_1, _mm256_load_ps(detail_avx2_dft8::orig_mult),
          _mm256_mul_ps(_mm256_permute_ps(temp_d_1, 0xB1),
                        _mm256_load_ps(detail_avx2_dft8::flip_mult)));

      _mm256_store_ps(out, _mm256_add_ps(d_0, d_1));
      _mm256_store_ps(out + 8, _mm256_sub_ps(d_0, d_1));
      MCA_END;
    }
  };
