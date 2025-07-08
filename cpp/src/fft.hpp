#include <array>
#include <complex>
#include <cstddef>
#include <cmath>
#include <arm_neon.h>

#define MCA_START __asm volatile("# LLVM-MCA-BEGIN");
#define MCA_END __asm volatile("# LLVM-MCA-END");

namespace FFT {
    // TODO : implement constexpr trigonometric functions

    // Neon intrinsics
    inline float32x2_t neon_mul(const float32x2_t &a, const float32x2_t &b) {
        return {a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
    }
    inline float neon_abs(const float32x2_t &a) {
        return std::sqrt(a[0]*a[0] + a[1]*a[1]);
    }
    inline float32x2_t neon_conj(const float32x2_t &a) {
        return {a[0], -a[1]};
    }
    
    constexpr double TwoPI = 2.0 * M_PI;
    constexpr double NegTwoPI = -TwoPI;
    

    // Just get the first smallest prime factor of N
    //  NOTE : this is key step in the induction for the mixed-radix-Cooley-Tukey FFT algorithm
    constexpr std::size_t get_prime_factor(std::size_t N) {
        for (std::size_t i = 2; i * i <= N; i++) {
            if (N % i == 0) {
                return i;
            }
        }
        return N;
    }

    // This is only for testing the get_prime_factor function
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

    // For a prime p of N, transpose the array so indexes of the same mod p are in order
    //  in the same "bin"
    //  NOTE : this is key step in the induction for the in-place mixed-radix-Cooley-Tukey FFT algorithm
    template <std::size_t N, typename T>
    constexpr void prime_factor_binner(T *data) {
        // Get first prime factor
        constexpr auto p = get_prime_factor(N);

        // Base case
        if constexpr (p < 2) {
            return;
        }

        // Have an array to track which elements have already been placed
        std::array<bool, N> flip_tracker = {false};

        // Given an index n, return the index it should be moved to
        constexpr auto get_next_index = [=](std::size_t n) -> std::size_t {
            const auto M = N / p;
            const auto i = n % p;
            const auto j = n / p;
            return i * M + j;
        };

        // Because of basic group theory, if N / p is not prime, we will
        //  have at least one loop of transposition that does not span all items.
        //  so as we whack-a-mole replace elements, we may need to find a new
        //  loop of elements to start replacing.
        std::size_t earliest_unmoved_index = 0;
        std::size_t curr_index = 0;
        std::size_t next_index = get_next_index(curr_index);
        T temp_val = data[curr_index];
        while (earliest_unmoved_index < N) {
            // Complete one loop of transpositions
            while (!flip_tracker[next_index]) {
                flip_tracker[next_index] = true;
                if (curr_index != next_index) { // avoid swap with self
                    std::swap(temp_val, data[next_index]);
                    curr_index = next_index;
                    next_index = get_next_index(curr_index);
                }
            }

            // Update earliest unmoved index if necessary
            while (earliest_unmoved_index < N && flip_tracker[earliest_unmoved_index])
                earliest_unmoved_index++;

            // Update bookkeeping for next loop of transpositions
            curr_index = earliest_unmoved_index;
            next_index = get_next_index(curr_index);
            temp_val = data[curr_index];
        }
        return;
    }

    // In-place Mixed Radix Cooley Tukey FFT
    template <std::size_t N>
    constexpr void fft_recurse(std::complex<float> *data) {
        // Get first prime factor
        constexpr auto p = get_prime_factor(N);
        constexpr auto M = N / p;

        // Base case
        if constexpr (p < 2) {
            return;
        }

        // Transpose Input Data around the radix
        MCA_START
        prime_factor_binner<N>(data);
        MCA_END

        // Recursively apply fft to each bin
        for (std::size_t i = 0; i < p; i++)
            fft_recurse<M>(data + (i * M));

        // Combine recursive results
        std::array<std::complex<float>, p> temp_factors;
        std::array<std::complex<float>, p> temp_results;
        constexpr auto get_factor = [](std::size_t i, std::size_t j, std::size_t k) -> std::complex<float> {
            const float angle = static_cast<float>(NegTwoPI * static_cast<double>(i * (k + j * M)) / static_cast<double>(N));
            return {std::cosf(angle), std::sinf(angle)};
        };
            
        for (std::size_t i = 0; i < M; i++) {
            // Store what strided values are
            for (std::size_t j = 0; j < p; j++)
                temp_factors[j] = data[i + (j * M)];

            // Recursively calculate strided results in temp location
            for (std::size_t j = 0; j < p; j++) {
                temp_results[j] = {0, 0};
                for (std::size_t k = 0; k < p; k++) {
                    temp_results[j] += temp_factors[k] * get_factor(k, j, i);
                }
            }
            
            // Finally write in strided results
            for (std::size_t j = 0; j < p; j++)
                data[i + (j * M)] = temp_results[j];
        }

        // This function is in-place! the input is modified with the result.
        return;
    }
    
    // In-place Mixed Radix Cooley Tukey FFT
    template <std::size_t N>
    constexpr void fft_recurse(float32x2_t *data) {
        // Get first prime factor
        constexpr auto p = get_prime_factor(N);
        constexpr auto M = N / p;

        // Base case
        if constexpr (p < 2) {
            return;
        }

        // Transpose Input Data around the radix
        prime_factor_binner<N>(data);

        // Recursively apply fft to each bin
        for (std::size_t i = 0; i < p; i++)
            fft_recurse<M>(data + (i * M));

        // Combine recursive results
        std::array<float32x2_t, p> temp_factors;
        std::array<float32x2_t, p> temp_results;
        
        // Setup constant factors for incremental rotations by multiplication
        constexpr float small_angle = static_cast<float>(NegTwoPI / static_cast<double>(N));
        constexpr float big_angle = static_cast<float>(NegTwoPI / static_cast<double>(p));
        const float32x2_t small_inc = {std::cosf(small_angle), std::sinf(small_angle)};
        const float32x2_t big_inc = {std::cosf(big_angle), std::sinf(big_angle)};

        float32x2_t i_factor = {1, 0};
        for (std::size_t i = 0; i < M; i++) {
            // Store what strided values are
            for (std::size_t j = 0; j < p; j++)
                temp_factors[j] = data[i + (j * M)];
        
            // Recursively calculate strided results in temp location
            float32x2_t j_factor = i_factor;
            for (std::size_t j = 0; j < p; j++) {
                temp_results[j] = {0, 0};
                float32x2_t k_factor = {1, 0};
                for (std::size_t k = 0; k < p; k++) {
                    // Finally, do the multiply and accumulate
                    temp_results[j] += neon_mul(temp_factors[k], k_factor);

                    // Increment k_factor
                    k_factor = neon_mul(k_factor, j_factor);
                }

                // Increment j_factor
                j_factor = neon_mul(j_factor, big_inc);
            }
            
            // Finally write in strided results
            for (std::size_t j = 0; j < p; j++)
                data[i + (j * M)] = temp_results[j];

            // Increment i_factor
            i_factor = neon_mul(i_factor, small_inc);
        }

        // This function is in-place! the input is modified with the result.
        return;
    }

    // TODO : may want to shift the result so that DC component is at the center,
    //  but without this shift the FFT and IFFT are inverses of each other
    template <std::size_t N>
    constexpr void fft(std::complex<float> *data) {
        fft_recurse<N>(data);

        // Normalize
        for (std::size_t i = 0; i < N; i++)
            data[i] /= static_cast<float>(N);

        return;
    }

    // TODO : may want to shift the result so that DC component is at the center,
    //  but without this shift the FFT and IFFT are inverses of each other
    template <std::size_t N>
    constexpr void fft(float32x2_t *data) {
        fft_recurse<N>(data);

        // Normalize
        for (std::size_t i = 0; i < N; i++)
            data[i] /= static_cast<float>(N);

        return;
    }

    template <std::size_t N>
    // Inverse FFT is just FFT on conj of input, and then conj of output
    constexpr void ifft(std::complex<float> *data) {
        for (std::size_t i = 0; i < N; i++)
            data[i] = std::conj(data[i]);

        fft_recurse<N>(data);

        for (std::size_t i = 0; i < N; i++)
            data[i] = std::conj(data[i]);

        // NOTE : if the fft is normalized, the ifft does not need to be normalized
        //  to keep the ifft the exact inverse operation of the fft.

        return;
    }
    
    template <std::size_t N>
    // Inverse FFT is just FFT on conj of input, and then conj of output
    constexpr void ifft(float32x2_t *data) {
        for (std::size_t i = 0; i < N; i++)
            data[i] = neon_conj(data[i]);

        fft_recurse<N>(data);

        for (std::size_t i = 0; i < N; i++)
            data[i] = neon_conj(data[i]);

        // NOTE : if the fft is normalized, the ifft does not need to be normalized
        //  to keep the ifft the exact inverse operation of the fft.

        return;
    }

    constexpr void wave_gen(std::complex<float> *time_domain, std::complex<float> *freq_domain,
        std::size_t N,
        unsigned int f = 1,
        unsigned int phase = 0,
        unsigned int amp = 1) {
        for (unsigned int i = 0; i < N; i++) {
            const float angle = static_cast<float>(TwoPI) * static_cast<float>((i * f) + phase)
                                  / static_cast<float>(N);
            time_domain[i] += std::complex<float>{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
        }
        const float angle = static_cast<float>(TwoPI) * 
                                  (static_cast<float>(phase) / static_cast<float>(N));
        freq_domain[f] += std::complex<float>{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
    }

    constexpr void wave_gen(float32x2_t *time_domain, float32x2_t *freq_domain,
        std::size_t N,
        unsigned int f = 1,
        unsigned int phase = 0,
        unsigned int amp = 1) {
        for (unsigned int i = 0; i < N; i++) {
            const float angle = static_cast<float>(TwoPI) * static_cast<float>((i * f) + phase)
                                  / static_cast<float>(N);
            time_domain[i] += float32x2_t{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
        }
        const float angle = static_cast<float>(TwoPI) * 
                                  (static_cast<float>(phase) / static_cast<float>(N));
        freq_domain[f] += float32x2_t{std::cosf(angle), std::sinf(angle)} * static_cast<float>(amp);
    }
}
