#include "fft.hpp"
// #include <fftw3.h>
#include <gperftools/profiler.h>

// place one of the benchmarks here to profile it
int main() {
    std::array<float32x2_t, 8192> time_domain1 = {0};
    std::array<float32x2_t, 8192> freq_domain1 = {0};
    for (std::size_t i = 0; i < 8192; i++) {
        FFT::wave_gen(time_domain1.data(), freq_domain1.data(), 8192,
            (i * 7) % 8191, (i + 3) % 7, i % 11);
    }
    ProfilerStart("fft power of 2");
    FFT::fft<8192>(time_domain1.data());
    ProfilerStop();
}