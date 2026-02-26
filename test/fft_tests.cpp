#include "../src/fft.hpp"
#include "../src/complex_types.hpp"
#include <fftw3.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>
#include <algorithm>

bool factors_check(unsigned int num, unsigned int *expected_factors, unsigned int expected_factor_count) {
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

// Your test cases go here
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

TEST_CASE("FFT Basic Small Input") {
    constexpr size_t N = 8;
    constexpr size_t M = 9;

    // Test Inputs
    std::array<std::complex<float>, N> time_domain1 = {0};
    std::array<std::complex<float>, N> freq_domain1 = {0};
    std::array<std::complex<float>, N> resp;
    std::array<std::complex<float>, N> resp2;
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    std::array<std::complex<float>, N> time_domain1_copy;
    for (std::size_t i = 0; i < N; i++)
        time_domain1_copy[i] = time_domain1[i];

    std::array<std::complex<float>, M> time_domain2 = {0};
    std::array<std::complex<float>, M> freq_domain2 = {0};
    FFT::wave_gen_lcg(time_domain2.data(), freq_domain2.data(), M);
    std::array<std::complex<float>, M> time_domain2_copy;
    for (std::size_t i = 0; i < M; i++)
        time_domain2_copy[i] = time_domain2[i];

    FFT::FFTPlan<FFT::StdComplexWrap<float>>::Init<N>();
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::Init<M>();

    // modify in place to frequency domain
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::fft<N>(
        static_cast<FFT::StdComplexWrap<float>*>(time_domain1.data()),
        static_cast<FFT::StdComplexWrap<float>*>(resp.data()));
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::fft<M>(
        static_cast<FFT::StdComplexWrap<float>*>(time_domain2.data()),
        static_cast<FFT::StdComplexWrap<float>*>(resp2.data()));

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_diff = std::max(max_diff, std::abs(resp[i] - freq_domain1[i]));
    }
    CHECK(max_diff < 2e-5);
    
    max_diff = 0.0f;
    for (std::size_t i = 0; i < M; i++) {
        max_diff = std::max(max_diff, std::abs(resp2[i] - freq_domain2[i]));
    }
    CHECK(max_diff < 6e-6);

    // modify in place back to time domain
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::ifft<N>(
        static_cast<FFT::StdComplexWrap<float>*>(resp.data()),
        static_cast<FFT::StdComplexWrap<float>*>(time_domain1.data()));
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::ifft<M>(
        static_cast<FFT::StdComplexWrap<float>*>(resp2.data()),
        static_cast<FFT::StdComplexWrap<float>*>(time_domain2.data()));
    
    // Check Outputs
    max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_diff = std::max(max_diff, std::abs(time_domain1[i] - time_domain1_copy[i]));
    }
    CHECK(max_diff < 8e-6);
    
    max_diff = 0;
    for (std::size_t i = 0; i < M; i++) {
        max_diff = std::max(max_diff, std::abs(time_domain2[i] - time_domain2_copy[i]));
    }
    CHECK(max_diff < 5e-6);
}

// Error is measured as the maximum magnitude of the difference
//  between any pair of computed and expected values
template<unsigned int N>
std::pair<float, float> fft_opt_tester() {
    // Test Inputs
    alignas(MY_MAX_ALIGNMENT) FFT::Complex time_domain[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex freq_domain[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex resp[N] = {0};
    FFT::wave_gen_lcg(time_domain, freq_domain, N);
    alignas(MY_MAX_ALIGNMENT) FFT::Complex time_domain_copy[N] = {0};
    for (std::size_t i = 0; i < N; i++)
        time_domain_copy[i] = time_domain[i];
    
    // Init

    FFT::FFTPlan<FFT::StdComplexWrap<float>>::Init<N>();
    
    // modify in place to frequency domain
    FFT::FFTPlan<FFT::Complex>::fft<N>(time_domain, resp);

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        /*std::cout << std::fixed << std::setprecision(2) << time_domain1[i][0] << " + i" << time_domain1[i][1] << " - " <<
          freq_domain1[i][0] << " + i" << freq_domain1[i][1] << std::endl;*/
        max_diff = std::max(max_diff, FFT::abs(resp[i] - freq_domain[i]));
    }

    // modify in place back to time domain
    FFT::FFTPlan<FFT::Complex>::ifft<N>(resp, time_domain);
    
    // Check Outputs
    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        /*std::cout << std::fixed << std::setprecision(2) << time_domain1[i][0] << " + i" << time_domain1[i][1] << " - " <<
          freq_domain1[i][0] << " + i" << freq_domain1[i][1] << std::endl;*/
        max_inverse_diff = std::max(max_inverse_diff, FFT::abs(time_domain[i] - time_domain_copy[i]));
    }
    return {max_diff, max_inverse_diff};
}

template<unsigned int N>
std::pair<float, float> fftw_tester() {
    // Configure input with waves
    fftwf_plan p;
    alignas(MY_MAX_ALIGNMENT) FFT::Complex in[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex out[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex in_cpy[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex out_cpy[N] = {0};
    for (std::size_t i = 0; i < N; i++) {
        in[i] = {0, 0};
        out[i] = {0, 0};
    }
    
    FFT::wave_gen_lcg(in, out, N);

    for (unsigned int i = 0; i < N; i++) {
        in_cpy[i] = in[i];
        out_cpy[i] = out[i];
    }

    // Set the plan
    p = fftwf_plan_dft_1d(N, 
        reinterpret_cast<fftwf_complex*>(in),
        reinterpret_cast<fftwf_complex*>(out),
        FFTW_FORWARD, FFTW_MEASURE);
    fftwf_execute(p);

    // Check Outputs
    float max_diff = 0;
    for (unsigned int i = 0; i < N; i++) {
        //std::cout << std::fixed << std::setprecision(2) << cpy[i].real() << " + i" << cpy[i].imag() << " - "
        //  << ans[i].real() << " + i" << ans[i].imag() << std::endl;
        max_diff = std::max(max_diff, FFT::abs(out[i] - static_cast<std::complex<float>>(N) * out_cpy[i]));
    }

    // Execute IFFT
    p = fftwf_plan_dft_1d(N, 
        reinterpret_cast<fftwf_complex*>(out),
        reinterpret_cast<fftwf_complex*>(in),
        FFTW_BACKWARD, FFTW_MEASURE);
    fftwf_execute(p);

    // Check Outputs
    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_inverse_diff = std::max(max_inverse_diff, FFT::abs(in[i] - static_cast<std::complex<float>>(N) * in_cpy[i]));
    }

    // Cleanup
    fftwf_destroy_plan(p);

    return {max_diff, max_inverse_diff};
}

// Try to ensure max error is ~1e-6
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
    CHECK(Ans_8.second < 2e-7);
}

// Keep error resonably small on small DFTs (un measurably small)
//  and not significant on medium DFTS (in the 1e-4 range)
TEST_CASE("FFT Opt Med Prime Input") {
    // Error on DFT of N 13 goes to float 0
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

// Tests FFTW compatability with reinterpret cast of std::complex<float>
//  and that it has similar overall error to custom FFT implementation
#ifdef ARCH_ARM
TEST_CASE("FFTW Compatability") {
    auto Ans_7 = fftw_tester<3 * 5 * 7>();
        CHECK(Ans_7.first < 7e-2);
        CHECK(Ans_7.second < 9e-3);
    }
#endif