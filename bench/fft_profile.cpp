#include "../src/complex_types.hpp"
#include "../src/fft.hpp"

// This is mostly used for perf and llvm-mca analysis
int main() {
    constexpr unsigned int N = 8192;
    alignas(MY_MAX_ALIGNMENT) FFT::Complex time_domain[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex freq_domain_ans[N] = {0};
    alignas(MY_MAX_ALIGNMENT) FFT::Complex freq_domain[N] = {0};
    FFT::wave_gen_lcg(time_domain, freq_domain_ans, N);

    FFT::FFTPlan<FFT::Complex>::Init<N>();
    for (unsigned int i = 0; i < 100000; i++) {
        FFT::FFTPlan<FFT::Complex>::fft<N>(time_domain, freq_domain);
    }
    return 0;
}