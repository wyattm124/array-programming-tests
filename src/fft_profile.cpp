#include "fft.hpp"
// #include <fftw3.h>
//#include <gperftools/profiler.h>

// place one of the benchmarks here to profile it
int main() {
    constexpr std::size_t N = 17;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);

    // TODO : have to initialize before Profile run!

    volatile auto fft_plan_n = FFT::FFTPlan<N>();
    //ProfilerStart("fft_mid_prime");
    for (std::size_t i = 0; i < 100000000; i++) {
        FFT::FFTPlan<N>::fft(time_domain1.data());
    }
    //ProfilerStop();
}