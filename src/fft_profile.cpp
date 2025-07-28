#include "fft.hpp"
// #include <fftw3.h>
//#include <gperftools/profiler.h>

// place one of the benchmarks here to profile it
int main() {
    constexpr std::size_t N = 13 * 17 * 7;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    for (std::size_t i = 0; i < N; i++) {
        FFT::wave_gen(time_domain1.data(), freq_domain1.data(), N,
            (i * 7) % (N-1)/* 8191 */, (i + 3) % 7, i % 11);
    }

    // TODO : have to initialize before Profile run!

    //ProfilerStart("fft power of 2");
    volatile auto fft_plan_n = FFT::FFTPlan<N>();
    FFT::FFTPlan<N>::fft(time_domain1.data());
    //ProfilerStop();
}