#include <iostream>
#include <iomanip>

constexpr unsigned int N = 7;

int main() {
    for (int i = 0; i < N; i++) {
        std::cout << std::fixed << std::setprecision(8) << std::cos((2 * M_PI * i) / N) 
        << " " << std::sin((2 * M_PI * i) / N) << std::endl;
    }    
}