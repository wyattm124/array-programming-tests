#include "cheb_approx.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

// Your test cases go here
 TEST_CASE("Cheb Poly") {
    bool cheb_poly_check = true;
    // T_0 = 1
    cheb_poly_check &= cheb_poly(0, 0) == 1;
    cheb_poly_check &= cheb_poly(0, 1) == 1;

    // T_1 = x
    cheb_poly_check &= cheb_poly(1, 0) == 0;
    cheb_poly_check &= cheb_poly(1, 1) == 1;
    cheb_poly_check &= cheb_poly(1, 2) == 2;

    // T_2 = 2x^2 - 1
    cheb_poly_check &= cheb_poly(2, -1) == 1;
    cheb_poly_check &= cheb_poly(2, 0) == -1;
    cheb_poly_check &= cheb_poly(2, 1) == 1;
    cheb_poly_check &= cheb_poly(2, 2) == 7;
    cheb_poly_check &= cheb_poly(2, 3) == 17;

    // t_3 = 4x^3 - 3x
    cheb_poly_check &= cheb_poly(3, -1) == -1;
    cheb_poly_check &= cheb_poly(3, 0) == 0;
    cheb_poly_check &= cheb_poly(3, 1) == 1;
    cheb_poly_check &= cheb_poly(3, 2) == 26;
    cheb_poly_check &= cheb_poly(3, 3) == 99;
    CHECK(cheb_poly_check);
 }

 TEST_CASE("Cheb Approx") {
    constexpr int range = 16;
    const auto poly1 = [](double x) -> double { return x; };
    const auto coef1 = cheb_gen_coef<2>(poly1);
    std::array<double, 2> coef1_array;
    std::copy(coef1.begin(), coef1.end(), coef1_array.begin());
    
    const auto poly2 = [](double x) -> double { return (2 * x * x) + (3 * x) + 1; };
    const auto coef2 = cheb_gen_coef<3>(poly2);
    std::array<double, 3> coef2_array;
    std::copy(coef2.begin(), coef2.end(), coef2_array.begin());
    
    const auto poly3 = [](double x) -> double { return (x * x * x) + (x * x) + (25 * x) + 17; };
    const auto coef3 = cheb_gen_coef<4>(poly3);
    std::array<double, 4> coef3_array;
    std::copy(coef3.begin(), coef3.end(), coef3_array.begin());

    double acc_error1 = 0;
    double acc_error2 = 0;
    double acc_error3 = 0;
    for (int i = -range; i < range; i++) {
        const double x = static_cast<double>(i) / range;
        acc_error1 += std::abs(poly1(x) - cheb_approx(x, coef1_array));
        acc_error2 += std::abs(poly2(x) - cheb_approx(x, coef2_array));
        acc_error3 += std::abs(poly3(x) - cheb_approx(x, coef3_array));
    }
    //CHECK(acc_error1 < 1e-6);
    //CHECK(acc_error2 < 1e-6);
    //CHECK(acc_error3 < 1e-6);
 }