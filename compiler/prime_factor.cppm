export module fft.compiler.prime_factor;

import <utility>;
import <vector>;

export namespace fft_compiler::prime_factor {

constexpr unsigned int get_prime_factor(unsigned int n) {
  for (unsigned int i = 2; i * i <= n; ++i) {
    if (n % i == 0) {
      return i;
    }
  }
  return n;
}

inline std::vector<std::pair<unsigned int, unsigned int>>
get_prime_factor_powers(unsigned int n) {
  std::vector<std::pair<unsigned int, unsigned int>> prime_factor_powers;

  if (n == 0) {
    return prime_factor_powers;
  }

  if (n == 1) {
    prime_factor_powers.emplace_back(1, 1);
    return prime_factor_powers;
  }

  while (n > 1) {
    const unsigned int prime_factor = get_prime_factor(n);
    unsigned int prime_exponent = 0;

    while (n % prime_factor == 0) {
      ++prime_exponent;
      n /= prime_factor;
    }

    prime_factor_powers.emplace_back(prime_factor, prime_exponent);
  }

  return prime_factor_powers;
}

} // namespace fft_compiler::prime_factor
