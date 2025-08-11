#include "fft.hpp"
#include <fftw3.h>
#include <benchmark/benchmark.h>


class Powerof2Fixture : public benchmark::Fixture {
protected:
    constexpr static size_t N = 8192;
    std::array<std::complex<float>, N> time_domain1;
    std::array<std::complex<float>, N> freq_domain1;
    
    void SetUp(const benchmark::State& state) {
        time_domain1.fill(0);
        freq_domain1.fill(0);
        FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    }
};

class MersennePrimeFixture : public benchmark::Fixture {
protected:
    constexpr static size_t N = 8191;
    std::array<std::complex<float>, N> time_domain1;
    std::array<std::complex<float>, N> freq_domain1;
        
    void SetUp(const benchmark::State& state) {
        time_domain1.fill(0);
        freq_domain1.fill(0);
        FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    }
};

class SmallPrimeFixture : public benchmark::Fixture {
protected:
    constexpr static size_t N = 17;
    std::array<std::complex<float>, N> time_domain1;
    std::array<std::complex<float>, N> freq_domain1;
        
    void SetUp(const benchmark::State& state) {
        time_domain1.fill(0);
        freq_domain1.fill(0);
        FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    }
};

// Power of 2 benchmarks
BENCHMARK_DEFINE_F(Powerof2Fixture, Basic)(benchmark::State& state) {
    for (auto _ : state) {
        FFT::FFTPlanBasic<Powerof2Fixture::N>::fft(time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(Powerof2Fixture, Neon)(benchmark::State& state) {
    // Need to Copy Data to Correct Type
    std::array<float32x2_t, Powerof2Fixture::N> temp_time_domain1 = {0};
    std::array<float32x2_t, Powerof2Fixture::N> temp_freq_domain1 = {0};
    for (unsigned int i = 0; i < Powerof2Fixture::N; i++) {
        temp_time_domain1[i] = {time_domain1[i].real(), time_domain1[i].imag()};
        temp_freq_domain1[i] = {freq_domain1[i].real(), freq_domain1[i].imag()};
    }

    // Make sure to do initialization
    volatile auto fft_plan_n = FFT::FFTPlan<Powerof2Fixture::N>();

    for (auto _ : state) {
        FFT::FFTPlan<Powerof2Fixture::N>::fft(temp_time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(Powerof2Fixture, FFTW)(benchmark::State& state) {
    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Powerof2Fixture::N);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Powerof2Fixture::N);
    for (std::size_t i = 0; i < Powerof2Fixture::N; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(Powerof2Fixture::N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

// Mersenne Prime benchmarks
BENCHMARK_DEFINE_F(MersennePrimeFixture, Basic)(benchmark::State& state) {
    for (auto _ : state) {
        FFT::FFTPlanBasic<MersennePrimeFixture::N>::fft(time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(MersennePrimeFixture, Neon)(benchmark::State& state) {
    // Need to Copy Data to Correct Type
    std::array<float32x2_t, MersennePrimeFixture::N> temp_time_domain1 = {0};
    std::array<float32x2_t, MersennePrimeFixture::N> temp_freq_domain1 = {0};
    for (unsigned int i = 0; i < MersennePrimeFixture::N; i++) {
        temp_time_domain1[i] = {time_domain1[i].real(), time_domain1[i].imag()};
        temp_freq_domain1[i] = {freq_domain1[i].real(), freq_domain1[i].imag()};
    }

    // Make sure to do initialization
    volatile auto fft_plan_n = FFT::FFTPlan<MersennePrimeFixture::N>();

    for (auto _ : state) {
        FFT::FFTPlan<MersennePrimeFixture::N>::fft(temp_time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(MersennePrimeFixture, FFTW)(benchmark::State& state) {
    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * MersennePrimeFixture::N);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * MersennePrimeFixture::N);
    for (std::size_t i = 0; i < MersennePrimeFixture::N; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(MersennePrimeFixture::N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

// Small Prime benchmarks
BENCHMARK_DEFINE_F(SmallPrimeFixture, Basic)(benchmark::State& state) {
    for (auto _ : state) {
        FFT::FFTPlanBasic<SmallPrimeFixture::N>::fft(time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(SmallPrimeFixture, Neon)(benchmark::State& state) {
    // Need to Copy Data to Correct Type
    std::array<float32x2_t, SmallPrimeFixture::N> temp_time_domain1 = {0};
    std::array<float32x2_t, SmallPrimeFixture::N> temp_freq_domain1 = {0};
    for (unsigned int i = 0; i < SmallPrimeFixture::N; i++) {
        temp_time_domain1[i] = {time_domain1[i].real(), time_domain1[i].imag()};
        temp_freq_domain1[i] = {freq_domain1[i].real(), freq_domain1[i].imag()};
    }

    // Make sure to do initialization
    volatile auto fft_plan_n = FFT::FFTPlan<SmallPrimeFixture::N>();

    for (auto _ : state) {
        FFT::FFTPlan<SmallPrimeFixture::N>::fft(temp_time_domain1.data());
    }
}

BENCHMARK_DEFINE_F(SmallPrimeFixture, FFTW)(benchmark::State& state) {
    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * SmallPrimeFixture::N);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * SmallPrimeFixture::N);
    for (std::size_t i = 0; i < SmallPrimeFixture::N; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(SmallPrimeFixture::N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

BENCHMARK_DEFINE_F(SmallPrimeFixture, Raders)(benchmark::State& state) {
    std::array<float32x2_t, 11> in;
    std::array<float32x2_t, 11> out = {0};
    for (unsigned int i = 0; i < 11; i++) {
        in[i] = {time_domain1[i].real(), time_domain1[i].imag()};
    }
    float32x2_t raders_coefs[7] = {
        {1.0f, 0.0f},
        {std::cosf((FFT::NegTwoPI / 7)), std::sinf((FFT::NegTwoPI / 7))},
        {std::cosf((2 * FFT::NegTwoPI / 7)), std::sinf((2 * FFT::NegTwoPI / 7))},
        {std::cosf((3 * FFT::NegTwoPI / 7)), std::sinf((3 * FFT::NegTwoPI / 7))},
        {std::cosf((4 * FFT::NegTwoPI / 7)), std::sinf((4 * FFT::NegTwoPI / 7))},
        {std::cosf((5 * FFT::NegTwoPI / 7)), std::sinf((5 * FFT::NegTwoPI / 7))},
        {std::cosf((6 * FFT::NegTwoPI / 7)), std::sinf((6 * FFT::NegTwoPI / 7))},
    };
    for (auto _ : state) {
        FFT::raders_7(in.data(), out.data(), raders_coefs);
    }
}

// Power of 2 Benchmarks
BENCHMARK_REGISTER_F(Powerof2Fixture, Basic);
BENCHMARK_REGISTER_F(Powerof2Fixture, Neon);
BENCHMARK_REGISTER_F(Powerof2Fixture, FFTW);

// Mersenne Prime Benchmarks
BENCHMARK_REGISTER_F(MersennePrimeFixture, Basic);
BENCHMARK_REGISTER_F(MersennePrimeFixture, Neon);
BENCHMARK_REGISTER_F(MersennePrimeFixture, FFTW);

// Small Prime Benchmarks
BENCHMARK_REGISTER_F(SmallPrimeFixture, Basic);
BENCHMARK_REGISTER_F(SmallPrimeFixture, Neon);
BENCHMARK_REGISTER_F(SmallPrimeFixture, FFTW);
BENCHMARK_REGISTER_F(SmallPrimeFixture, Raders);

// Template out main
BENCHMARK_MAIN();