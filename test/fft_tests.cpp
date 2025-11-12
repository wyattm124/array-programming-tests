#include "../src/fft.hpp"
#include "../src/complex_types.hpp"
#include <fftw3.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>

// Your test cases go here
TEST_CASE("Prime Factorization") {
    CHECK(prime_factor::prime_factorization(1) == std::array<unsigned int, 64>{1});
    CHECK(prime_factor::prime_factorization(2) == std::array<unsigned int, 64>{2});
    CHECK(prime_factor::prime_factorization(3) == std::array<unsigned int, 64>{3});
    CHECK(prime_factor::prime_factorization(4) == std::array<unsigned int, 64>{2, 2});
    CHECK(prime_factor::prime_factorization(5) == std::array<unsigned int, 64>{5});
    CHECK(prime_factor::prime_factorization(6) == std::array<unsigned int, 64>{2, 3});
    CHECK(prime_factor::prime_factorization(7) == std::array<unsigned int, 64>{7});
    CHECK(prime_factor::prime_factorization(8) == std::array<unsigned int, 64>{2, 2, 2});
    CHECK(prime_factor::prime_factorization(9) == std::array<unsigned int, 64>{3, 3});
    CHECK(prime_factor::prime_factorization(10) == std::array<unsigned int, 64>{2, 5});
    CHECK(prime_factor::prime_factorization(11) == std::array<unsigned int, 64>{11});
    CHECK(prime_factor::prime_factorization(12) == std::array<unsigned int, 64>{2, 2, 3});
    CHECK(prime_factor::prime_factorization(13) == std::array<unsigned int, 64>{13});
    CHECK(prime_factor::prime_factorization(14) == std::array<unsigned int, 64>{2, 7});
    CHECK(prime_factor::prime_factorization(15) == std::array<unsigned int, 64>{3, 5});
    CHECK(prime_factor::prime_factorization(16) == std::array<unsigned int, 64>{2, 2, 2, 2});
    CHECK(prime_factor::prime_factorization(17) == std::array<unsigned int, 64>{17});
    CHECK(prime_factor::prime_factorization(18) == std::array<unsigned int, 64>{2, 3, 3});
    CHECK(prime_factor::prime_factorization(19) == std::array<unsigned int, 64>{19});
    CHECK(prime_factor::prime_factorization(20) == std::array<unsigned int, 64>{2, 2, 5});
    CHECK(prime_factor::prime_factorization(21) == std::array<unsigned int, 64>{3, 7});
    CHECK(prime_factor::prime_factorization(22) == std::array<unsigned int, 64>{2, 11});
    CHECK(prime_factor::prime_factorization(23) == std::array<unsigned int, 64>{23});
    CHECK(prime_factor::prime_factorization(24) == std::array<unsigned int, 64>{2, 2, 2, 3});
    CHECK(prime_factor::prime_factorization(25) == std::array<unsigned int, 64>{5, 5});
    CHECK(prime_factor::prime_factorization(26) == std::array<unsigned int, 64>{2, 13});
    CHECK(prime_factor::prime_factorization(27) == std::array<unsigned int, 64>{3, 3, 3});
    CHECK(prime_factor::prime_factorization(28) == std::array<unsigned int, 64>{2, 2, 7});
    CHECK(prime_factor::prime_factorization(29) == std::array<unsigned int, 64>{29});
    CHECK(prime_factor::prime_factorization(30) == std::array<unsigned int, 64>{2, 3, 5});
    CHECK(prime_factor::prime_factorization(10007) == std::array<unsigned int, 64>{10007});
    CHECK(prime_factor::prime_factorization(10008) == std::array<unsigned int, 64>{2, 2, 2, 3, 3, 139});
    CHECK(prime_factor::prime_factorization(13 * 17) == std::array<unsigned int, 64>{13, 17});
    CHECK(prime_factor::prime_factorization(13 * 17 * 19) == std::array<unsigned int, 64>{13, 17, 19});
}

TEST_CASE("FFT Radix Transposition Extra Memory") {
    // Test Inputs
    std::array<float32x2_t, 8> data1 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}}};
    std::array<float32x2_t, 9> data2 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}, {8, 0}}};
    std::array<float32x2_t, 8> out1;
    std::array<float32x2_t, 9> out2;

    // Test Answers
    std::array<float32x2_t, 8> answer1 = {{{0, 0}, {2, 0}, {4, 0}, {6, 0}, {1, 0}, {3, 0}, {5, 0}, {7, 0}}};
    std::array<float32x2_t, 9> answer2 = {{{0, 0}, {3, 0}, {6, 0}, {1, 0}, {4, 0}, {7, 0}, {2, 0}, {5, 0}, {8, 0}}};
    
    // modify in place
    FFT::DFT_binner<4, 2>(data1.data(), out1.data());
    FFT::DFT_binner<3, 3>(data2.data(), out2.data());
    
    // Check Outputs
    bool cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (out1[i][0] == answer1[i][0] && out1[i][1] == answer1[i][1]);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (out2[i][0] == answer2[i][0] && out2[i][1] == answer2[i][1]);
    }
    CHECK(cmp);
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

    // modify in place to frequency domain
    FFT::FFTPlan<N, FFT::StdComplexWrap<float>>::fft(
        static_cast<FFT::StdComplexWrap<float>*>(time_domain1.data()),
        static_cast<FFT::StdComplexWrap<float>*>(resp.data()));
    FFT::FFTPlan<M, FFT::StdComplexWrap<float>>::fft(
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
    FFT::FFTPlan<N, FFT::StdComplexWrap<float>>::ifft(
        static_cast<FFT::StdComplexWrap<float>*>(resp.data()),
        static_cast<FFT::StdComplexWrap<float>*>(time_domain1.data()));
    FFT::FFTPlan<M, FFT::StdComplexWrap<float>>::ifft(
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
    std::array<FFT::Complex, N> time_domain = {0};
    std::array<FFT::Complex, N> freq_domain = {0};
    std::array<FFT::Complex, N> resp;
    FFT::wave_gen_lcg(time_domain.data(), freq_domain.data(), N);
    std::array<FFT::Complex, N> time_domain_copy;
    for (std::size_t i = 0; i < N; i++)
        time_domain_copy[i] = time_domain[i];
    
    // Init
    volatile auto fft_plan = FFT::FFTPlan<N, FFT::Complex>();
    
    // modify in place to frequency domain
    FFT::FFTPlan<N, FFT::Complex>::fft(time_domain.data(), resp.data());

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        /*std::cout << std::fixed << std::setprecision(2) << time_domain1[i][0] << " + i" << time_domain1[i][1] << " - " <<
          freq_domain1[i][0] << " + i" << freq_domain1[i][1] << std::endl;*/
        max_diff = std::max(max_diff, FFT::abs(resp[i] - freq_domain[i]));
    }

    // modify in place back to time domain
    FFT::FFTPlan<N, FFT::Complex>::ifft(resp.data(), time_domain.data());
    
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
    std::array<std::complex<float>, N> in;
    std::array<std::complex<float>, N> out;
    std::array<std::complex<float>, N> in_cpy;
    std::array<std::complex<float>, N> out_cpy;
    for (std::size_t i = 0; i < N; i++) {
        in[i] = {0, 0};
        out[i] = {0, 0};
    }
    
    FFT::wave_gen_lcg(in.data(), out.data(), N);

    for (std::size_t i = 0; i < N; i++) {
        in_cpy[i] = in[i];
        out_cpy[i] = out[i];
    }

    // Set the plan
    p = fftwf_plan_dft_1d(N, 
        reinterpret_cast<fftwf_complex*>(in.data()),
        reinterpret_cast<fftwf_complex*>(out.data()),
        FFTW_FORWARD, FFTW_MEASURE);
    fftwf_execute(p);

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        //std::cout << std::fixed << std::setprecision(2) << cpy[i].real() << " + i" << cpy[i].imag() << " - "
        //  << ans[i].real() << " + i" << ans[i].imag() << std::endl;
        max_diff = std::max(max_diff, std::abs(out[i] - static_cast<std::complex<float>>(N) * out_cpy[i]));
    }

    // Execute IFFT
    p = fftwf_plan_dft_1d(N, 
        reinterpret_cast<fftwf_complex*>(out.data()),
        reinterpret_cast<fftwf_complex*>(in.data()),
        FFTW_BACKWARD, FFTW_MEASURE);
    fftwf_execute(p);

    // Check Outputs
    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_inverse_diff = std::max(max_inverse_diff, std::abs(in[i] - static_cast<std::complex<float>>(N) * in_cpy[i]));
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
    CHECK(Ans_5.second < 3e-7);
    
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
    CHECK(Ans_53.second < 6.2e-5);
}

TEST_CASE("FFT Opt Small Prime Composite Input") {
    auto Ans_1 = fft_opt_tester<3 * 5>();
    CHECK(Ans_1.first < 2e-5);
    CHECK(Ans_1.second < 2e-2);
    
    auto Ans_2 = fft_opt_tester<7 * 5>();
    CHECK(Ans_2.first < 9e-5);
    CHECK(Ans_2.second < 3e-5);

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
TEST_CASE("FFTW Compatability") {
    auto Ans_7 = fftw_tester<3 * 5 * 7>();
    CHECK(Ans_7.first < 7e-2);
    CHECK(Ans_7.second < 9e-3);
}