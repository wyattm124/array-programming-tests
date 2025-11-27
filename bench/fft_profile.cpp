#include "../src/fft.hpp"
#include "fft_comp_unit.hpp"
//#include <fftw3.h>
//#include <gperftools/profiler.h>

// This is mostly used to see the instructions used to implement an FFT case
int main() {
    constexpr std::size_t N = 5;
    std::array<FFT::Complex, N> time_domain = {0};
    std::array<FFT::Complex, N> freq_domain_ans = {0};
    std::array<FFT::Complex, N> freq_domain = {0};
    FFT::wave_gen_lcg(time_domain.data(), freq_domain_ans.data(), N);

    FFT::FFTPlan<FFT::Complex>::Init<N>();
    //ProfilerStart("fft_mid_prime");
    for (std::size_t i = 0; i < 100000000; i++) {
        FFT::FFTPlan<FFT::Complex>::fft<N>(time_domain.data(), freq_domain.data());
    }
    //ProfilerStop();
}