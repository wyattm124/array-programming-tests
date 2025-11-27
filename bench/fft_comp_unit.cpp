#include "../src/fft.hpp"
#include "fft_comp_unit.hpp"

namespace FFT {
    void fft_2(Complex* __restrict__ in, Complex* __restrict__ out) noexcept { 
        FFT::FFTPlan<Complex>::fft<2>(in, out);
    }
    void fft_3(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<3>(in, out);
    }
    void fft_4(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<4>(in, out);
    }
    void fft_5(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<5>(in, out);
    }
    void fft_6(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<6>(in, out);
    }
    void fft_7(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<7>(in, out);
    }
    void fft_8(Complex* __restrict__ in, Complex* __restrict__ out) noexcept {
        FFT::FFTPlan<Complex>::fft<8>(in, out);
    }
}
        