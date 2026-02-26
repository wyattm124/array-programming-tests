#include <iostream>
#include <iomanip>
#include <cmath>

constexpr unsigned int N = 8;

// Does not add accuracy for Complex float calculations to 
//   approximate roots past 8 digits of precision
int main() {
    for (int i = 0; i < N; i++) {
        std::cout << std::fixed << std::setprecision(12) << cos((2 * M_PI * i) / N) 
        << " " << sin((2 * M_PI * i) / N) << std::endl;
    }    
}