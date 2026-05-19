#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "../src/fft.hpp"
#include <algorithm>
#include <array>
#include <fftw3.h>
#include <utility>

bool factors_check(unsigned int num, unsigned int *expected_factors,
                   unsigned int expected_factor_count) {
    unsigned int factors[64];
    prime_factor::prime_factorization(num, factors);
    for (unsigned int i = 0; i < 64; i++) {
        if (i < expected_factor_count) {
            if (factors[i] != expected_factors[i])
                return false;
        } else {
            if (factors[i] != 0)
                return false;
        }
    }
    return true;
}

template <typename T>
T pair_abs(T real, T imag) {
    return static_cast<T>(std::sqrt(real * real + imag * imag));
}

template <typename T, std::size_t N>
T max_pair_diff(const std::array<T, N> &lhs, const std::array<T, N> &rhs) {
    T max_diff = 0;
    for (std::size_t i = 0; i < N / 2; i++) {
        max_diff = std::max(
            max_diff,
            pair_abs(lhs[2 * i] - rhs[2 * i], lhs[2 * i + 1] - rhs[2 * i + 1]));
    }
    return max_diff;
}

template <unsigned int N>
void interleaved_to_split(const float (&interleaved)[2 * N], float (&real)[N],
                          float (&imag)[N]) {
    for (unsigned int i = 0; i < N; i++) {
        real[i] = interleaved[2 * i];
        imag[i] = interleaved[2 * i + 1];
    }
}

template <unsigned int N>
void split_to_interleaved(const float (&real)[N], const float (&imag)[N],
                          float (&interleaved)[2 * N]) {
    for (unsigned int i = 0; i < N; i++) {
        interleaved[2 * i] = real[i];
        interleaved[2 * i + 1] = imag[i];
    }
}

TEST_CASE("Prime Factorization") {
    unsigned int ans_1[1] = {1};
    CHECK(factors_check(1, ans_1, 1));

    unsigned int ans_2[1] = {2};
    CHECK(factors_check(2, ans_2, 1));

    unsigned int ans_5[1] = {5};
    CHECK(factors_check(5, ans_5, 1));

    unsigned int ans_6[2] = {2, 3};
    CHECK(factors_check(6, ans_6, 2));

    unsigned int ans_12[3] = {2, 2, 3};
    CHECK(factors_check(12, ans_12, 3));

    unsigned int ans_16[4] = {2, 2, 2, 2};
    CHECK(factors_check(16, ans_16, 4));

    unsigned int ans_17[1] = {17};
    CHECK(factors_check(17, ans_17, 1));

    unsigned int ans_24[4] = {2, 2, 2, 3};
    CHECK(factors_check(24, ans_24, 4));

    unsigned int ans_26[2] = {2, 13};
    CHECK(factors_check(26, ans_26, 2));

    unsigned int ans_28[3] = {2, 2, 7};
    CHECK(factors_check(28, ans_28, 3));

    unsigned int ans_10007[1] = {10007};
    CHECK(factors_check(10007, ans_10007, 1));

    unsigned int ans_10008[6] = {2, 2, 2, 3, 3, 139};
    CHECK(factors_check(10008, ans_10008, 6));

    unsigned int ans_13_17[2] = {13, 17};
    CHECK(factors_check(13 * 17, ans_13_17, 2));

    unsigned int ans_13_17_19[3] = {13, 17, 19};
    CHECK(factors_check(13 * 17 * 19, ans_13_17_19, 3));
}

template<unsigned int N, SIMD_TYPE simd = SIMD_TYPE::NONE>
std::pair<float, float> fft_opt_tester() {
    alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float resp[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float time_domain_copy[2 * N] = {0};
    FFT::wave_gen_lcg(time_domain, freq_domain, N);
    for (std::size_t i = 0; i < 2 * N; i++)
        time_domain_copy[i] = time_domain[i];

    FFT::FFTPlan<float, simd>::template Init<N>();
    FFT::FFTPlan<float, simd>::template fft<N>(time_domain, resp);

    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_diff = std::max(max_diff,
            pair_abs(resp[2 * i] - freq_domain[2 * i],
                     resp[2 * i + 1] - freq_domain[2 * i + 1]));
    }

    FFT::FFTPlan<float, simd>::template ifft<N>(resp, time_domain);

    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_inverse_diff = std::max(max_inverse_diff,
            pair_abs(time_domain[2 * i] - time_domain_copy[2 * i],
                     time_domain[2 * i + 1] - time_domain_copy[2 * i + 1]));
    }
    return {max_diff, max_inverse_diff};
}

template<unsigned int N>
std::pair<float, float> fftw_tester() {
    fftwf_plan p;
    fftwf_iodim dims[1] = {{static_cast<int>(N), 1, 1}};
    alignas(MY_MAX_ALIGNMENT) float in[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float out[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float in_cpy[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float out_cpy[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float in_real[N] = {0};
    alignas(MY_MAX_ALIGNMENT) float in_imag[N] = {0};
    alignas(MY_MAX_ALIGNMENT) float out_real[N] = {0};
    alignas(MY_MAX_ALIGNMENT) float out_imag[N] = {0};

    FFT::wave_gen_lcg(in, out, N);

    for (unsigned int i = 0; i < 2 * N; i++) {
        in_cpy[i] = in[i];
        out_cpy[i] = out[i];
    }

    interleaved_to_split(in, in_real, in_imag);
    p = fftwf_plan_guru_split_dft(1, dims, 0, nullptr,
        in_real, in_imag, out_real, out_imag, FFTW_MEASURE);
    fftwf_execute(p);
    split_to_interleaved(out_real, out_imag, out);
    fftwf_destroy_plan(p);

    float max_diff = 0;
    for (unsigned int i = 0; i < N; i++) {
        max_diff = std::max(max_diff,
            pair_abs(out[2 * i] - static_cast<float>(N) * out_cpy[2 * i],
                     out[2 * i + 1] - static_cast<float>(N) * out_cpy[2 * i + 1]));
    }

    interleaved_to_split(out, out_real, out_imag);
    p = fftwf_plan_guru_split_dft(1, dims, 0, nullptr,
        out_real, out_imag, in_real, in_imag, FFTW_MEASURE);
    fftwf_execute(p);
    split_to_interleaved(in_real, in_imag, in);

    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_inverse_diff = std::max(max_inverse_diff,
            pair_abs(in[2 * i] - static_cast<float>(N) * in_cpy[2 * i],
                     in[2 * i + 1] - static_cast<float>(N) * in_cpy[2 * i + 1]));
    }
    fftwf_destroy_plan(p);

    return {max_diff, max_inverse_diff};
}

TEST_CASE("FFT Basic Small Input") {
    auto Ans_8 = fft_opt_tester<8>();
    CHECK(Ans_8.first < 2e-5f);
    CHECK(Ans_8.second < 8e-6f);

    auto Ans_9 = fft_opt_tester<9>();
    CHECK(Ans_9.first < 6e-6f);
    CHECK(Ans_9.second < 5e-6f);
}

TEST_CASE("FFT Opt Base Case Input") {
    auto Ans_3 = fft_opt_tester<3>();
    CHECK(Ans_3.first < 2e-7);
    CHECK(Ans_3.second < 2e-7);

    auto Ans_4 = fft_opt_tester<4>();
    CHECK(Ans_4.first < 5e-7);
    CHECK(Ans_4.second < 2e-7);

    auto Ans_5 = fft_opt_tester<5>();
    CHECK(Ans_5.first < 6e-7);
    CHECK(Ans_5.second < 3.5e-7);

    auto Ans_6 = fft_opt_tester<6>();
    CHECK(Ans_6.first < 1.3e-6);
    CHECK(Ans_6.second < 3e-7);

    auto Ans_7 = fft_opt_tester<7>();
    CHECK(Ans_7.first < 9.3e-7);
    CHECK(Ans_7.second < 4e-7);

    auto Ans_8 = fft_opt_tester<8>();
    CHECK(Ans_8.first < 1.1e-6);
    CHECK(Ans_8.second < 5e-7);
}

TEST_CASE("FFT Opt Med Prime Input") {
    auto Ans_13 = fft_opt_tester<13>();
    CHECK(Ans_13.first < 1e-9);
    CHECK(Ans_13.second < 1e-9);

    auto Ans_53 = fft_opt_tester<53>();
    CHECK(Ans_53.first < 2.3e-4);
    CHECK(Ans_53.second < 7e-5);
}

TEST_CASE("FFT Opt Small Prime Composite Input") {
    auto Ans_1 = fft_opt_tester<3 * 5>();
    CHECK(Ans_1.first < 2e-5);
    CHECK(Ans_1.second < 2e-2);

    auto Ans_2 = fft_opt_tester<7 * 5>();
    CHECK(Ans_2.first < 9e-5);
    CHECK(Ans_2.second < 3.5e-5);

    auto Ans_3 = fft_opt_tester<3 * 5 * 5>();
    CHECK(Ans_3.first < 4e-4);
    CHECK(Ans_3.second < 3e-1);

    auto Ans_4 = fft_opt_tester<3 * 5 * 7>();
    CHECK(Ans_4.first < 7e-4);
    CHECK(Ans_4.second < 3e-1);

    auto Ans_5 = fft_opt_tester<7 * 11>();
    CHECK(Ans_5.first < 4e-4);
    CHECK(Ans_5.second < 6e-5);

    auto Ans_6 = fft_opt_tester<3 * 5 * 7 * 11 * 13>();
    CHECK(Ans_6.first < 2e-3);
    CHECK(Ans_6.second < 2e-2);
}

#ifdef ARCH_ARM
TEST_CASE("FFTW Compatability") {
    auto Ans_7 = fftw_tester<3 * 5 * 7>();
    CHECK(Ans_7.first < 7e-2);
    CHECK(Ans_7.second < 9e-3);
}
#endif
