#include "../src/prime_factor.hpp"
#include <cstddef>
#include <iostream>

namespace prime_factor {
    template <std::size_t N>
    bool n_smooth(std::size_t num) {
        constexpr std::size_t primes[64] = {
            2, 3, 5, 7, 11, 13, 17, 19,
            23, 29, 31, 37, 41, 43, 47, 53,
            59, 61, 67, 71, 73, 79, 83, 89,
            97, 101, 103, 107, 109, 113, 127, 131,
            137, 139, 149, 151, 157, 163, 167, 173,
            179, 181, 191, 193, 197, 199, 211, 223,
            227, 229, 233, 239, 241, 251, 257, 263,
            269, 271, 277, 281, 283, 293, 307, 311
        };
        for (std::size_t i = 0; i < N; i++) {
            while (num % primes[i] == 0) { num /= primes[i]; }
        }
        return num == 1;
    }
    template <std::size_t N>
    std::size_t n_smooth_dist(std::size_t num) {
        std::size_t last_n_smooth = 1;
        std::size_t max_dist = 0;
        for (std::size_t i = 2; i <= num; i++) {
            if (n_smooth<N>(i)) {
                max_dist = std::max(max_dist, i - last_n_smooth);
                last_n_smooth = i;
            }
        }
        return max_dist;
    }
}

// for 2^32 takes about 1 min - note, this is all the way to about 4.2 billion
// max 29 smooth dist is 423,752
// max 31 smooth dist is 201,175
// max 37 smooth dist is 182,584
int main() {
    std::cout << prime_factor::n_smooth_dist<12>(4294967295) << std::endl;
    return 0;
}
    