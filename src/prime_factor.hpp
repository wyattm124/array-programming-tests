#pragma once

namespace prime_factor {
    // Just get the first smallest prime factor of N
    // NOTE : this is key step in the induction for the mixed-radix-Cooley-Tukey FFT algorithm
    constexpr std::size_t get_prime_factor(std::size_t N) {
        for (std::size_t i = 2; i * i <= N; i++) {
            if (N % i == 0) {
                return i;
            }
        }
        return N;
    }

    // NOTE : This is only for testing the get_prime_factor function
    constexpr std::array<std::size_t, 64> prime_factorization(std::size_t n) {
        // Take care of "degenerate" edge cases
        std::array<std::size_t, 64> factors = {0};
        if (n == 0) {
            return factors;
        } else if (n == 1) {
            factors[0] = 1;
            return factors;
        }

        // For all other cases, do the usual loop
        std::size_t factor_count = 0;
        std::size_t p = get_prime_factor(n);
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