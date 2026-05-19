#include "../src/fft.hpp"
#include "fft_comp_unit.hpp"

namespace FFT {
    void fft_2(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<2>(in, out);
    }
    void fft_3(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<3>(in, out);
    }
    void fft_4(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<4>(in, out);
    }
    void fft_5(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<5>(in, out);
    }
    void fft_6(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<6>(in, out);
    }
    void fft_7(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<7>(in, out);
    }
    void fft_8(float* __restrict__ in, float* __restrict__ out) noexcept {
        FFT::FFTPlan<float>::fft<8>(in, out);
    }
}
