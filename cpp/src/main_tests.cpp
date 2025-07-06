/* Restrictions:
 * (1) all data accessed by these operations must be of known size of compile time known
 *   size and type
 * (2) all indexing operations must be constexpr qualifiable - that is they can be fully
 *   fully defined and potentially evaluted at compile time
 * (3) 
 */

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

 TEST_CASE("FFT Radix Transposition") {
    // Test Inputs
    std::array<std::complex<float>, 8> data1 = {{0, 1, 2, 3, 4, 5, 6, 7}};
    std::array<std::complex<float>, 9> data2 = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};

    // Test Answers
    std::array<std::complex<float>, 8> answer1 = {{0, 2, 4, 6, 1, 3, 5, 7}};
    std::array<std::complex<float>, 9> answer2 = {{0, 3, 6, 1, 4, 7, 2, 5, 8}};
    
    // modify in place
    FFT::prime_factor_binner<8>(data1.data());
    FFT::prime_factor_binner<9>(data2.data());
    
    // Check Outputs
    for (std::size_t i = 0; i < 8; i++) {
        CHECK(data1[i] == answer1[i]);
    }
    for (std::size_t i = 0; i < 9; i++) {
        CHECK(data2[i] == answer2[i]);
    }
 }

 TEST_CASE("FFT") {
    // Test Inputs
    std::array<std::complex<float>, 8> data1 = {0};
    FFT::wave_gen(data1.data(), 8, 1, 0, 1);
    FFT::wave_gen(data1.data(), 8, 2, 0, 1);

    std::array<std::complex<float>, 9> data2 = {0};
    FFT::wave_gen(data2.data(), 9, 1, 0, 1);
    FFT::wave_gen(data2.data(), 9, 3, 0, 1);

    // TODO : Test Answers

    // modify in place
    FFT::fft<8>(data1.data());
    FFT::fft<9>(data2.data());
    
    // Check Outputs
    std::cout << "Case 1" << std::endl;
    for (std::size_t i = 0; i < 8; i++) {
        std::cout << data1[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << "Case 2" << std::endl;
    for (std::size_t i = 0; i < 9; i++) {
        std::cout << data2[i] << ", ";
    }
    std::cout << std::endl;
 }