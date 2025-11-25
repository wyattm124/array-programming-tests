#pragma once

#include <array>

namespace prime_factor {
    // Just get the first smallest prime factor of N
    // NOTE : this is key step in the induction for the mixed-radix-Cooley-Tukey FFT algorithm
    constexpr unsigned int get_prime_factor(unsigned int N) {
        for (unsigned int i = 2; i * i <= N; i++) {
            if (N % i == 0) {
                return i;
            }
        }
        return N;
    }

    template<unsigned int n>
    concept Prime = (get_prime_factor(n) == n);

    // NOTE : This is only for testing the get_prime_factor function
    constexpr std::array<unsigned int, 64> prime_factorization(unsigned int n) {
        // Take care of "degenerate" edge cases
        std::array<unsigned int, 64> factors = {0};
        if (n == 0) {
            return factors;
        } else if (n == 1) {
            factors[0] = 1;
            return factors;
        }

        // For all other cases, do the usual loop
        unsigned int factor_count = 0;
        unsigned int p = get_prime_factor(n);
        while (p > 1) {
            while (n % p == 0) {
                factors[factor_count++] = p;
                n /= p;
            }
            p = get_prime_factor(n);
        }
        return factors;
    }
}