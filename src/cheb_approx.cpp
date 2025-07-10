#include <iostream>
#include <vector>
#include <functional>

/*
 * This was made based on Wolfram Mathworld's article on Chebyshev Approximation
 * https://mathworld.wolfram.com/ChebyshevApproximationFormula.html
 * 
 * std::cos and std::sin are not constexpr functions in the llvm toolchain used
 * so we need to use constexpr versions of these functions. Polynomials can be
 * evaluated at compile time using constexpr, so the best way to approximate
 * std::sin and std::cos are to use Chebyshev polynomials approximations.
 * 
 * As std::cos and std:;sin need to be approximated over a large range, a Maclaurin
 * series or a Taylor series polynomial at any point is not a good option as they
 * do not converge well over a large range.
 *
 * As the functions to be approximated need evaluated over a long range, lookup tables
 * are not practical and will not perform well either.
 *
 */

constexpr std::size_t factorial(std::size_t n) noexcept {
    if (n == 0) return 1;
    return n * factorial(n - 1);
}

constexpr std::size_t combinitorial(std::size_t n, std::size_t k) noexcept {
    if (k > n || k == 0 || n == 0) return 1;
    return factorial(n) / (factorial(k) * factorial(n - k));
}

constexpr double cheb_poly(double in, std::size_t N) noexcept {
    // factors to grow each term by
    const double a = in * in;
    const double b = a - 1;

    // accumulating coefficients
    double acc_a = 1;
    double acc_b = 1;
    for (std::size_t i = 0; i < N; i++)
        acc_a *= a;

    // Compute the chebyshev polynomial as a series
    double result = 0;
    for (std::size_t i = 0; i < N/2; i++) {
        result += combinitorial(N, 2 * i) * acc_a * acc_b;
        acc_a /= a;
        acc_b *= b;
    }
    return result;
}

std::vector<double> cheb_gen_coef(std::size_t N, const std::function<double(double)> &f) noexcept {
    std::vector<double> result;
    for (std::size_t i = 0; i < N; i++) {
      double coef = 0;
      for (std::size_t j = 1; j < N + 1; j++) {
        const double x = std::cos(M_PI * (static_cast<double>(i) - 0.5) / static_cast<double>(N));
        coef += f(x) * cheb_poly(x, i);
      }
      result.push_back((2 * coef) / static_cast<double>(N));
    }
    return result;
}

template <std::size_t N>
double cheb_approx(double in, const std::array<double, N> &coef) noexcept {
    double result = -0.5 * coef[0];
    for (std::size_t i = 0; i < N; i++) {
        result += coef[i] * cheb_poly(in, i);
    }
    return result;
}

// Put any expression in this main, and it can be approximated using Chebyshev polynomials
//  with the coefficients being written to a CSV.
int main() {
    // TODO : write tests to make sure this approximates well
    // TODO : generate coefficients for sin and cos
    std::cout << "Hello, World!" << std::endl;
    return 0;
}