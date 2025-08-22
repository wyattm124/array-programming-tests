#include <iostream>
#include <iomanip>

constexpr unsigned int N = 8;

// does not add accuracy past 8 digits of precision
int main() {
    for (int i = 0; i < N; i++) {
        std::cout << std::fixed << std::setprecision(12) << std::cos((2 * M_PI * i) / N) 
        << " " << std::sin((2 * M_PI * i) / N) << std::endl;
    }    
}