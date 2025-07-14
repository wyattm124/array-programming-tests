#include "fft.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

// Your test cases go here
TEST_CASE("Prime Factorization") {
    CHECK(FFT::prime_factorization(1) == std::array<size_t, 64>{1});
    CHECK(FFT::prime_factorization(2) == std::array<size_t, 64>{2});
    CHECK(FFT::prime_factorization(3) == std::array<size_t, 64>{3});
    CHECK(FFT::prime_factorization(4) == std::array<size_t, 64>{2, 2});
    CHECK(FFT::prime_factorization(5) == std::array<size_t, 64>{5});
    CHECK(FFT::prime_factorization(6) == std::array<size_t, 64>{2, 3});
    CHECK(FFT::prime_factorization(7) == std::array<size_t, 64>{7});
    CHECK(FFT::prime_factorization(8) == std::array<size_t, 64>{2, 2, 2});
    CHECK(FFT::prime_factorization(9) == std::array<size_t, 64>{3, 3});
    CHECK(FFT::prime_factorization(10) == std::array<size_t, 64>{2, 5});
    CHECK(FFT::prime_factorization(11) == std::array<size_t, 64>{11});
    CHECK(FFT::prime_factorization(12) == std::array<size_t, 64>{2, 2, 3});
    CHECK(FFT::prime_factorization(13) == std::array<size_t, 64>{13});
    CHECK(FFT::prime_factorization(14) == std::array<size_t, 64>{2, 7});
    CHECK(FFT::prime_factorization(15) == std::array<size_t, 64>{3, 5});
    CHECK(FFT::prime_factorization(16) == std::array<size_t, 64>{2, 2, 2, 2});
    CHECK(FFT::prime_factorization(17) == std::array<size_t, 64>{17});
    CHECK(FFT::prime_factorization(18) == std::array<size_t, 64>{2, 3, 3});
    CHECK(FFT::prime_factorization(19) == std::array<size_t, 64>{19});
    CHECK(FFT::prime_factorization(20) == std::array<size_t, 64>{2, 2, 5});
    CHECK(FFT::prime_factorization(21) == std::array<size_t, 64>{3, 7});
    CHECK(FFT::prime_factorization(22) == std::array<size_t, 64>{2, 11});
    CHECK(FFT::prime_factorization(23) == std::array<size_t, 64>{23});
    CHECK(FFT::prime_factorization(24) == std::array<size_t, 64>{2, 2, 2, 3});
    CHECK(FFT::prime_factorization(25) == std::array<size_t, 64>{5, 5});
    CHECK(FFT::prime_factorization(26) == std::array<size_t, 64>{2, 13});
    CHECK(FFT::prime_factorization(27) == std::array<size_t, 64>{3, 3, 3});
    CHECK(FFT::prime_factorization(28) == std::array<size_t, 64>{2, 2, 7});
    CHECK(FFT::prime_factorization(29) == std::array<size_t, 64>{29});
    CHECK(FFT::prime_factorization(30) == std::array<size_t, 64>{2, 3, 5});
    CHECK(FFT::prime_factorization(10007) == std::array<size_t, 64>{10007});
    CHECK(FFT::prime_factorization(10008) == std::array<size_t, 64>{2, 2, 2, 3, 3, 139});
    CHECK(FFT::prime_factorization(13 * 17) == std::array<size_t, 64>{13, 17});
    CHECK(FFT::prime_factorization(13 * 17 * 19) == std::array<size_t, 64>{13, 17, 19});
}

TEST_CASE("FFT Radix Transposition Opt Inplace") {
    // Test Inputs
    std::array<float32x2_t, 8> in1 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}}};
    std::array<float32x2_t, 9> in2 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}, {8, 0}}};

    // Test Answers
    std::array<float32x2_t, 8> answer1 = {{{0, 0}, {2, 0}, {4, 0}, {6, 0}, {1, 0}, {3, 0}, {5, 0}, {7, 0}}};
    std::array<float32x2_t, 9> answer2 = {{{0, 0}, {3, 0}, {6, 0}, {1, 0}, {4, 0}, {7, 0}, {2, 0}, {5, 0}, {8, 0}}};
    
    // modify in place
    FFT::FFTPlan<8>::prime_factor_binner_inplace(in1.data());
    FFT::FFTPlan<9>::prime_factor_binner_inplace(in2.data());
    
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

TEST_CASE("FFT Radix Transposition Opt Extra Memory") {
    // Test Inputs
    std::array<float32x2_t, 8> data1 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}}};
    std::array<float32x2_t, 9> data2 = {{{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}, {8, 0}}};
    std::array<float32x2_t, 8> out1;
    std::array<float32x2_t, 9> out2;

    // Test Answers
    std::array<float32x2_t, 8> answer1 = {{{0, 0}, {2, 0}, {4, 0}, {6, 0}, {1, 0}, {3, 0}, {5, 0}, {7, 0}}};
    std::array<float32x2_t, 9> answer2 = {{{0, 0}, {3, 0}, {6, 0}, {1, 0}, {4, 0}, {7, 0}, {2, 0}, {5, 0}, {8, 0}}};
    
    // modify in place
    FFT::FFTPlan<8>::prime_factor_binner(data1.data(), out1.data());
    FFT::FFTPlan<9>::prime_factor_binner(data2.data(), out2.data());
    
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

TEST_CASE("FFT Radix Transposition Basic") {
    // Test Inputs
    std::array<std::complex<float>, 8> data1 = {{0, 1, 2, 3, 4, 5, 6, 7}};
    std::array<std::complex<float>, 9> data2 = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};
    std::array<std::complex<float>, 8> out1;
    std::array<std::complex<float>, 9> out2;

    // Test Answers
    std::array<std::complex<float>, 8> answer1 = {{0, 2, 4, 6, 1, 3, 5, 7}};
    std::array<std::complex<float>, 9> answer2 = {{0, 3, 6, 1, 4, 7, 2, 5, 8}};
    
    // modify in place
    FFT::FFTPlanBasic<8>::prime_factor_binner(data1.data(), out1.data());
    FFT::FFTPlanBasic<9>::prime_factor_binner(data2.data(), out2.data());
    
    // Check Outputs
    bool cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (out1[i] == answer1[i]);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (out2[i] == answer2[i]);
    }
    CHECK(cmp);
}

TEST_CASE("FFT Basic Small Input") {
    // Test Inputs
    std::array<std::complex<float>, 8> time_domain1 = {0};
    std::array<std::complex<float>, 8> freq_domain1 = {0};
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), 8, 1, 1, 1);
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), 8, 2, 3, 1);
    std::array<std::complex<float>, 8> time_domain1_copy;
    for (std::size_t i = 0; i < 8; i++)
        time_domain1_copy[i] = time_domain1[i];

    std::array<std::complex<float>, 9> time_domain2 = {0};
    std::array<std::complex<float>, 9> freq_domain2 = {0};
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), 9, 1, 3, 1);
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), 9, 3, 4, 1);
    std::array<std::complex<float>, 9> time_domain2_copy;
    for (std::size_t i = 0; i < 9; i++)
        time_domain2_copy[i] = time_domain2[i]; 

    // modify in place to frequency domain
    FFT::FFTPlanBasic<8>::fft(time_domain1.data());
    FFT::FFTPlanBasic<9>::fft(time_domain2.data());

    // Check Outputs
    bool cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (std::abs(time_domain1[i] - freq_domain1[i]) < 1e-6);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (std::abs(time_domain2[i] - freq_domain2[i]) < 1e-6);
    }
    CHECK(cmp);

    // modify in place back to time domain
    FFT::FFTPlanBasic<8>::ifft(time_domain1.data());
    FFT::FFTPlanBasic<9>::ifft(time_domain2.data());
    
    // Check Outputs
    cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (std::abs(time_domain1[i] - time_domain1_copy[i]) < 1e-6);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (std::abs(time_domain2[i] - time_domain2_copy[i]) < 1e-6);
    }
    CHECK(cmp);
}

 TEST_CASE("FFT Opt Small Input") {
    // Test Inputs
    std::array<float32x2_t, 8> time_domain1 = {0};
    std::array<float32x2_t, 8> freq_domain1 = {0};
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), 8, 1, 1, 1);
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), 8, 2, 3, 1);
    std::array<float32x2_t, 8> time_domain1_copy;
    for (std::size_t i = 0; i < 8; i++)
        time_domain1_copy[i] = time_domain1[i];

    std::array<float32x2_t, 9> time_domain2 = {0};
    std::array<float32x2_t, 9> freq_domain2 = {0};
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), 9, 1, 3, 1);
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), 9, 3, 4, 1);
    std::array<float32x2_t, 9> time_domain2_copy;
    for (std::size_t i = 0; i < 9; i++)
        time_domain2_copy[i] = time_domain2[i]; 

    // modify in place to frequency domain
    FFT::FFTPlan<8>::fft(time_domain1.data());
    FFT::FFTPlan<9>::fft(time_domain2.data());

    // Check Outputs
    bool cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (FFT::neon_abs(time_domain1[i] - freq_domain1[i]) < 1e-6);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (FFT::neon_abs(time_domain2[i] - freq_domain2[i]) < 1e-6);
    }
    CHECK(cmp);

    // modify in place back to time domain
    FFT::FFTPlan<8>::ifft(time_domain1.data());
    FFT::FFTPlan<9>::ifft(time_domain2.data());
    
    // Check Outputs
    cmp = true;
    for (std::size_t i = 0; i < 8; i++) {
        cmp &= (FFT::neon_abs(time_domain1[i] - time_domain1_copy[i]) < 1e-6);
    }
    CHECK(cmp);
    
    cmp = true;
    for (std::size_t i = 0; i < 9; i++) {
        cmp &= (FFT::neon_abs(time_domain2[i] - time_domain2_copy[i]) < 1e-6);
    }
    CHECK(cmp);
}

// TODO : these are outside bounds of error!
TEST_CASE("FFT Opt Med Input") {
    constexpr size_t N = 128;//8192;
    constexpr size_t M = 127;//8191;

    // Test Inputs
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), N, 1, 1, 1);
    FFT::wave_gen(time_domain1.data(), freq_domain1.data(), N, 2, 3, 1);
    std::array<float32x2_t, N> time_domain1_copy;
    for (std::size_t i = 0; i < N; i++)
        time_domain1_copy[i] = time_domain1[i];

    std::array<float32x2_t, M> time_domain2 = {0};
    std::array<float32x2_t, M> freq_domain2 = {0};
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), M, 1, 3, 1);
    FFT::wave_gen(time_domain2.data(), freq_domain2.data(), M, 3, 4, 1);
    std::array<float32x2_t, M> time_domain2_copy;
    for (std::size_t i = 0; i < M; i++)
        time_domain2_copy[i] = time_domain2[i]; 

    // modify in place to frequency domain
    FFT::FFTPlan<N>::fft(time_domain1.data());
    FFT::FFTPlan<M>::fft(time_domain2.data());

    // Check Outputs
    float max_error = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_error = std::max(max_error, FFT::neon_abs(time_domain1[i] - freq_domain1[i]));
    }
    CHECK(max_error < 5e-7);
    
    max_error = 0;
    for (std::size_t i = 0; i < M; i++) {
        max_error = std::max(max_error, FFT::neon_abs(time_domain2[i] - freq_domain2[i]));
    }
    CHECK(max_error < 5e-7);

    // modify in place back to time domain
    FFT::FFTPlan<N>::ifft(time_domain1.data());
    FFT::FFTPlan<M>::ifft(time_domain2.data());
    
    // Check Outputs
    max_error = 0;
    for (std::size_t i = 0; i < N; i++) {
        max_error = std::max(max_error, FFT::neon_abs(time_domain1[i] - time_domain1_copy[i]));
    }
    CHECK(max_error < 5e-7);
    
    max_error = 0;
    for (std::size_t i = 0; i < M; i++) {
        max_error = std::max(max_error, FFT::neon_abs(time_domain2[i] - time_domain2_copy[i]));
    }
    CHECK(max_error < 2e-6);
}