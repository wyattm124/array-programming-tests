#include "../src/fft.hpp"
#include "fft_comp_unit.hpp"

namespace FFT {
    void fft_2(Complex* __restrict__ in, Complex* __restrict__ out) noexcept { 
        FFT::FFTPlan<2, Complex>::fft(in, out);
    }
    void fft_3(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<3, Complex>::fft(in, out);
    }
    void fft_4(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<4, Complex>::fft(in, out);
    }
    void fft_5(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<5, Complex>::fft(in, out);
    }
    void fft_6(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<6, Complex>::fft(in, out);
    }
    void fft_7(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<7, Complex>::fft(in, out);
    }
    void fft_8(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<8, Complex>::fft(in, out);
    }
}
        