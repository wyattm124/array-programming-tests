#include "fft.hpp"
// #include <fftw3.h>
//#include <gperftools/profiler.h>

// place one of the benchmarks here to profile it
int main() {
    constexpr std::size_t N = 8191;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    for (std::size_t i = 0; i < N; i++) {
        FFT::wave_gen(time_domain1.data(), freq_domain1.data(), N,
            (i * 7) % 8191, (i + 3) % 7, i % 11);
    }
    //ProfilerStart("fft power of 2");
    FFT::fft<N>(time_domain1.data());
    //ProfilerStop();
}