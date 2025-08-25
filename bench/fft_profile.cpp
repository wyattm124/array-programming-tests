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

    volatile auto fft_plan_n = FFT::FFTPlan<N, FFT::Complex>();
    //ProfilerStart("fft_mid_prime");
    for (std::size_t i = 0; i < 100000000; i++) {
        FFT::FFTPlan<N, FFT::Complex>::fft(time_domain.data(), freq_domain.data());
    }
    //ProfilerStop();
}