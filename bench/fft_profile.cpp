#include "../src/fft.hpp"

int main() {
    constexpr unsigned int N = 8192;
    alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float freq_domain_ans[2 * N] = {0};
    alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N] = {0};
    FFT::wave_gen_lcg(time_domain, freq_domain_ans, N);

    FFT::FFTPlan<float>::Init<N>();
    for (unsigned int i = 0; i < 100000; i++) {
        FFT::FFTPlan<float>::fft<N>(time_domain, freq_domain);
    }
    return 0;
}
