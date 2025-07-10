#include <iostream>

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

// Put any expression in this main, and it can be approximated using Chebyshev polynomials
//  with the coefficients being written to a CSV.
int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}