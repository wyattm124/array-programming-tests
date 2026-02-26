#pragma once

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
    constexpr void prime_factorization(unsigned int n, unsigned int (&factors)[64]) {
        // Make sure to clear the factors array
        for (unsigned int i = 0; i < 64; i++) {
            factors[i] = 0;
        }

        // Take care of "degenerate" edge cases
        if (n == 0) {
            return;
        } else if (n == 1) {
            factors[0] = 1;
            return;
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
        return;
    }
}