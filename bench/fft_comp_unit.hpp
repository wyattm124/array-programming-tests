#pragma once

#include "../src/complex_types.hpp"

namespace FFT {
    void fft_2(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_3(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_4(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_5(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_6(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_7(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
    void fft_8(Complex* __restrict__ in, Complex* __restrict__ out) noexcept;
}