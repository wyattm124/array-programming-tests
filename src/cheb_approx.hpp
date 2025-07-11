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

// Use the recurrence relation to compute the Chebyshev polynomials,
//   and if possible do this at compile time
constexpr double cheb_poly(std::size_t N, double x) noexcept {
    if (N == 0) return 1;
    double n_2 = 1;
    double n_1 = x;
    double n_0 = x; // not quite correct, but helps case N = 1
    for (std::size_t i = 2; i <= N; i++) {
        n_0 = 2 * x * n_1 - n_2;
        n_2 = n_1;
        n_1 = n_0;
    }
    return n_0;
}

template <std::size_t N>
constexpr std::array<double, N> cheb_gen_coef(const std::function<double(double)> &f) noexcept {
    std::array<double, N> result;
    for (std::size_t i = 0; i < N; i++) {
      double coef = 0;
      for (std::size_t j = 1; j < N + 1; j++) {
        const double x = std::cos(M_PI * (static_cast<double>(i) - 0.5) / static_cast<double>(N));
        coef += f(x) * cheb_poly(i, x);
      }
      result[i] = (2 * coef) / static_cast<double>(N);
    }
    return result;
}

template <std::size_t N>
constexpr double cheb_approx(double x, const std::array<double, N> &coef) noexcept {
    double result = -0.5 * coef[0];
    for (std::size_t i = 1; i < N; i++) {
        result += coef[i] * cheb_poly(i, x);
    }
    return result;
}