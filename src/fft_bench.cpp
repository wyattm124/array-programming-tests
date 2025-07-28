#include "fft.hpp"
#include <fftw3.h>
#include <benchmark/benchmark.h>

// Define another benchmark
static void BM_Basic_Powerof2FFT(benchmark::State& state) {
    std::array<std::complex<float>, 8192> time_domain1 = {0};
    std::array<std::complex<float>, 8192> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), 8192);
    for (auto _ : state) {
        FFT::FFTPlanBasic<8192>::fft(time_domain1.data());
    }
}

// Define another benchmark
static void BM_Neon_Powerof2FFT(benchmark::State& state) {
    constexpr size_t N = 8192;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);

    // Make sure coeffs are calculated ahead of time
    volatile auto fft_plan_n = FFT::FFTPlan<N>();
    for (auto _ : state) {
        FFT::FFTPlan<N>::fft(time_domain1.data());
    }
}

static void BM_FFTW_Powerof2FFT(benchmark::State& state) {
    // Generate Waves
    std::array<std::complex<float>, 8192> time_domain1 = {0};
    std::array<std::complex<float>, 8192> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), 8192);

    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 8192);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 8192);
    for (std::size_t i = 0; i < 8192; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(8192, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

static void BM_Basic_MersennePrimeFFT(benchmark::State& state) {
    std::array<std::complex<float>, 8191> time_domain1 = {0};
    std::array<std::complex<float>, 8191> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), 8191);
    for (auto _ : state) {
        FFT::FFTPlanBasic<8191>::fft(time_domain1.data());
    }
}

static void BM_Neon_MersennePrimeFFT(benchmark::State& state) {
    constexpr size_t N = 8191;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);

    // Make sure coeffs are calculated ahead of time
    volatile auto fft_plan_n = FFT::FFTPlan<N>();
    for (auto _ : state) {
        FFT::FFTPlan<N>::fft(time_domain1.data());
    }
}

static void BM_FFTW_MersennePrimeFFT(benchmark::State& state) {
    // Generate Waves
    std::array<std::complex<float>, 8191> time_domain1 = {0};
    std::array<std::complex<float>, 8191> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), 8191);

    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 8191);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 8191);
    for (std::size_t i = 0; i < 8191; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(8191, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

static void BM_Basic_SmallPrimeFFT(benchmark::State& state) {
    constexpr size_t N = 17;
    std::array<std::complex<float>, N> time_domain1 = {0};
    std::array<std::complex<float>, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);
    for (auto _ : state) {
        FFT::FFTPlanBasic<N>::fft(time_domain1.data());
    }
}

static void BM_Neon_SmallPrimeFFT(benchmark::State& state) {
    constexpr size_t N = 17;
    std::array<float32x2_t, N> time_domain1 = {0};
    std::array<float32x2_t, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);

    // Make sure coeffs are calculated ahead of time
    volatile auto fft_plan_n = FFT::FFTPlan<N>();
    for (auto _ : state) {
        FFT::FFTPlan<N>::fft(time_domain1.data());
    }
}

static void BM_FFTW_SmallPrimeFFT(benchmark::State& state) {
    constexpr size_t N = 17;
    // Generate Waves
    std::array<std::complex<float>, N> time_domain1 = {0};
    std::array<std::complex<float>, N> freq_domain1 = {0};
    FFT::wave_gen_lcg(time_domain1.data(), freq_domain1.data(), N);

    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
    for (std::size_t i = 0; i < N; i++) {
        in[i][0] = time_domain1[i].real();
        in[i][1] = time_domain1[i].imag();
        out[i][0] = 0;
        out[i][1] = 0;
    }

    // Set the plan
    p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

// Add all Benchmarks
BENCHMARK(BM_Basic_Powerof2FFT);
BENCHMARK(BM_Neon_Powerof2FFT);
BENCHMARK(BM_FFTW_Powerof2FFT);
BENCHMARK(BM_Basic_MersennePrimeFFT);
BENCHMARK(BM_Neon_MersennePrimeFFT);
BENCHMARK(BM_FFTW_MersennePrimeFFT);
BENCHMARK(BM_Basic_SmallPrimeFFT);
BENCHMARK(BM_Neon_SmallPrimeFFT);
BENCHMARK(BM_FFTW_SmallPrimeFFT);

// Template out main
BENCHMARK_MAIN();