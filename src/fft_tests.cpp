#include "fft.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <iomanip>

// Your test cases go here
TEST_CASE("Prime Factorization") {
    CHECK(prime_factor::prime_factorization(1) == std::array<size_t, 64>{1});
    CHECK(prime_factor::prime_factorization(2) == std::array<size_t, 64>{2});
    CHECK(prime_factor::prime_factorization(3) == std::array<size_t, 64>{3});
    CHECK(prime_factor::prime_factorization(4) == std::array<size_t, 64>{2, 2});
    CHECK(prime_factor::prime_factorization(5) == std::array<size_t, 64>{5});
    CHECK(prime_factor::prime_factorization(6) == std::array<size_t, 64>{2, 3});
    CHECK(prime_factor::prime_factorization(7) == std::array<size_t, 64>{7});
    CHECK(prime_factor::prime_factorization(8) == std::array<size_t, 64>{2, 2, 2});
    CHECK(prime_factor::prime_factorization(9) == std::array<size_t, 64>{3, 3});
    CHECK(prime_factor::prime_factorization(10) == std::array<size_t, 64>{2, 5});
    CHECK(prime_factor::prime_factorization(11) == std::array<size_t, 64>{11});
    CHECK(prime_factor::prime_factorization(12) == std::array<size_t, 64>{2, 2, 3});
    CHECK(prime_factor::prime_factorization(13) == std::array<size_t, 64>{13});
    CHECK(prime_factor::prime_factorization(14) == std::array<size_t, 64>{2, 7});
    CHECK(prime_factor::prime_factorization(15) == std::array<size_t, 64>{3, 5});
    CHECK(prime_factor::prime_factorization(16) == std::array<size_t, 64>{2, 2, 2, 2});
    CHECK(prime_factor::prime_factorization(17) == std::array<size_t, 64>{17});
    CHECK(prime_factor::prime_factorization(18) == std::array<size_t, 64>{2, 3, 3});
    CHECK(prime_factor::prime_factorization(19) == std::array<size_t, 64>{19});
    CHECK(prime_factor::prime_factorization(20) == std::array<size_t, 64>{2, 2, 5});
    CHECK(prime_factor::prime_factorization(21) == std::array<size_t, 64>{3, 7});
    CHECK(prime_factor::prime_factorization(22) == std::array<size_t, 64>{2, 11});
    CHECK(prime_factor::prime_factorization(23) == std::array<size_t, 64>{23});
    CHECK(prime_factor::prime_factorization(24) == std::array<size_t, 64>{2, 2, 2, 3});
    CHECK(prime_factor::prime_factorization(25) == std::array<size_t, 64>{5, 5});
    CHECK(prime_factor::prime_factorization(26) == std::array<size_t, 64>{2, 13});
    CHECK(prime_factor::prime_factorization(27) == std::array<size_t, 64>{3, 3, 3});
    CHECK(prime_factor::prime_factorization(28) == std::array<size_t, 64>{2, 2, 7});
    CHECK(prime_factor::prime_factorization(29) == std::array<size_t, 64>{29});
    CHECK(prime_factor::prime_factorization(30) == std::array<size_t, 64>{2, 3, 5});
    CHECK(prime_factor::prime_factorization(10007) == std::array<size_t, 64>{10007});
    CHECK(prime_factor::prime_factorization(10008) == std::array<size_t, 64>{2, 2, 2, 3, 3, 139});
    CHECK(prime_factor::prime_factorization(13 * 17) == std::array<size_t, 64>{13, 17});
    CHECK(prime_factor::prime_factorization(13 * 17 * 19) == std::array<size_t, 64>{13, 17, 19});
}

TEST_CASE("FFT Radix Transposition Inplace") {
    // Test Inputs
    std::array<float32x2_t, 8> in1 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}}};
    std::array<float32x2_t, 9> in2 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}, {8, 0}}};

    // Test Answers
    std::array<float32x2_t, 8> answer1 = {{{0, 0}, {2, 0}, {4, 0}, {6, 0}, {1, 0}, {3, 0}, {5, 0}, {7, 0}}};
    std::array<float32x2_t, 9> answer2 = {{{0, 0}, {3, 0}, {6, 0}, {1, 0}, {4, 0}, {7, 0}, {2, 0}, {5, 0}, {8, 0}}};
    
    // modify in place
    FFT::DFT_binner_inplace<4, 2>(in1.data());
    FFT::DFT_binner_inplace<3, 3>(in2.data());
    
    // Check Outputs
    bool cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (in1[i][0] == answer1[i][0] && in1[i][1] == answer1[i][1]);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (in2[i][0] == answer2[i][0] && in2[i][1] == answer2[i][1]);
    }
    CHECK(cmp);
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
    FFT::FFTPlanBasic<N>::fft(time_domain1.data());
    FFT::FFTPlanBasic<M>::fft(time_domain2.data());

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_diff = std::max(max_diff, std::abs(time_domain1[i] - freq_domain1[i]));
    }
    CHECK(max_diff < 2e-5);
    
    max_diff = 0.0f;
    for (std::size_t i = 0; i < M; i++) {
        max_diff = std::max(max_diff, std::abs(time_domain2[i] - freq_domain2[i]));
    }
    CHECK(max_diff < 6e-6);

    // modify in place back to time domain
    FFT::FFTPlanBasic<N>::ifft(time_domain1.data());
    FFT::FFTPlanBasic<M>::ifft(time_domain2.data());
    
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

template<unsigned int N>
std::pair<float, float> fft_opt_tester() {
    // Test Inputs
    std::array<FFT::Complex, N> time_domain1 = {0};
    std::array<FFT::Complex, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    std::array<FFT::Complex, N> time_domain1_copy;
    for (std::size_t i = 0; i < N; i++)
        time_domain1_copy[i] = time_domain1[i];
    
    // Init
    volatile auto fft_plan = FFT::FFTPlan<N, FFT::Complex>();
    
    // modify in place to frequency domain
    FFT::FFTPlan<N, FFT::Complex>::fft(time_domain1.data());

    // Check Outputs
    float max_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        /*std::cout << std::fixed << std::setprecision(2) << time_domain1[i][0] << " + i" << time_domain1[i][1] << " - " <<
          freq_domain1[i][0] << " + i" << freq_domain1[i][1] << std::endl;*/
        max_diff = std::max(max_diff, FFT::abs(time_domain1[i] - freq_domain1[i]));
    }

    // modify in place back to time domain
    FFT::FFTPlan<N, FFT::Complex>::ifft(time_domain1.data());
    
    // Check Outputs
    float max_inverse_diff = 0;
    for (std::size_t i = 0; i < N; i++) {
        /*std::cout << std::fixed << std::setprecision(2) << time_domain1[i][0] << " + i" << time_domain1[i][1] << " - " <<
          freq_domain1[i][0] << " + i" << freq_domain1[i][1] << std::endl;*/
        max_inverse_diff = std::max(max_inverse_diff, FFT::abs(time_domain1[i] - time_domain1_copy[i]));
    }
    return {max_diff, max_inverse_diff};
}
TEST_CASE("FFT Opt Base Case Input") {
    auto Ans_3 = fft_opt_tester<3>();
    CHECK(Ans_3.first < 4e-4);
    CHECK(Ans_3.second < 2e-3);

    auto Ans_4 = fft_opt_tester<4>();
    CHECK(Ans_4.first < 2e-5);
    CHECK(Ans_4.second < 3e-6);
    
    auto Ans_5 = fft_opt_tester<5>();
    CHECK(Ans_5.first < 8e-6);
    CHECK(Ans_5.second < 3e-6);
    
    auto Ans_6 = fft_opt_tester<6>();
    CHECK(Ans_6.first < 8e-6);
    CHECK(Ans_6.second < 3e-6);

    auto Ans_7 = fft_opt_tester<7>();
    CHECK(Ans_7.first < 1e-6);
    CHECK(Ans_7.second < 4e-6);

    auto Ans_8 = fft_opt_tester<8>();
    CHECK(Ans_8.first < 2e-6);
    CHECK(Ans_8.second < 5e-7);
}

TEST_CASE("FFT Opt Med Prime Input") {
    auto Ans_13 = fft_opt_tester<13>();
    CHECK(Ans_13.first < 8e-6);
    CHECK(Ans_13.second < 2e-6);
    
    auto Ans_53 = fft_opt_tester<53>();
    CHECK(Ans_53.first < 3e-4);
    CHECK(Ans_53.second < 7e-5);
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

// TODO : Opt tests with larger numbers like 8192 and 8191?