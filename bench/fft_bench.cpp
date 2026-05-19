#include "fft_comp_unit.hpp"
#include "../src/fft.hpp"
#include <benchmark/benchmark.h>
#include <fftw3.h>

#include <algorithm>
#include <vector>

template <unsigned long M>
void interleaved_to_split(const float (&interleaved)[M], float *real,
                          float *imag) {
    static_assert((M % 2) == 0, "Interleaved complex buffers must have even length.");
    constexpr unsigned long N = M / 2;
    for (unsigned long i = 0; i < N; i++) {
        real[i] = interleaved[2 * i];
        imag[i] = interleaved[2 * i + 1];
    }
}

/// Base Class for Sizes
template <unsigned int M>
class Case : public benchmark::Fixture {
public:
    constexpr static unsigned int N = M;
protected:
    alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N];
    alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N];

    void SetUp(const benchmark::State& state) {
        (void)state;
        std::fill_n(time_domain, 2 * N, 0.0f);
        std::fill_n(freq_domain, 2 * N, 0.0f);
        FFT::wave_gen_lcg(time_domain, freq_domain, N);
    }
};

using Case2 = Case<2>;
using Case3 = Case<3>;
using Case4 = Case<4>;
using Case5 = Case<5>;
using Case6 = Case<6>;
using Case7 = Case<7>;
using Case8 = Case<8>;
using PowerOf2 = Case<8192>;
using MediumPrime = Case<53>;
using MersennePrime = Case<8191>;

template <unsigned long M>
void process_fftw(benchmark::State& state, float (&time_domain)[M]) {
    static_assert((M % 2) == 0, "Interleaved complex buffers must have even length.");
    constexpr unsigned long N = M / 2;
    fftwf_iodim dims[1] = {{static_cast<int>(N), 1, 1}};
    float *in_real = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
    float *in_imag = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
    float *out_real = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
    float *out_imag = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
    interleaved_to_split(time_domain, in_real, in_imag);
    fftwf_plan p = fftwf_plan_guru_split_dft(
        1, dims, 0, nullptr, in_real, in_imag, out_real, out_imag, FFTW_MEASURE);

    for (auto _ : state) {
        fftwf_execute(p);
    }

    fftwf_destroy_plan(p);
    fftwf_free(in_real);
    fftwf_free(in_imag);
    fftwf_free(out_real);
    fftwf_free(out_imag);
}

template <unsigned long M>
void process_fftw_with_alloc(benchmark::State& state, float (&time_domain)[M]) {
    static_assert((M % 2) == 0, "Interleaved complex buffers must have even length.");
    constexpr unsigned long N = M / 2;
    fftwf_iodim dims[1] = {{static_cast<int>(N), 1, 1}};
    fftwf_plan p;
    std::vector<float *> buffs;

    for (auto _ : state) {
        float *in_real = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
        float *in_imag = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
        float *out_real = static_cast<float*>(fftwf_malloc(sizeof(float) * N));
        float *out_imag = static_cast<float*>(fftwf_malloc(sizeof(float) * N));

        p = fftwf_plan_guru_split_dft(
            1, dims, 0, nullptr, in_real, in_imag, out_real, out_imag, FFTW_ESTIMATE);

        interleaved_to_split(time_domain, in_real, in_imag);
        fftwf_execute(p);

        buffs.push_back(in_real);
        buffs.push_back(in_imag);
        buffs.push_back(out_real);
        buffs.push_back(out_imag);
        fftwf_destroy_plan(p);
    }

    for (auto *buf : buffs) {
        fftwf_free(buf);
    }
}

template <unsigned long M>
void process_fft(benchmark::State& state, float (&time_domain)[M]) {
    static_assert((M % 2) == 0, "Interleaved complex buffers must have even length.");
    constexpr unsigned long N = M / 2;
    alignas(MY_MAX_ALIGNMENT) float temp_freq_domain[M] = {0};
    FFT::FFTPlan<float>::Init<N>();

    for (auto _ : state) {
        FFT::FFTPlan<float>::fft<N>(time_domain, temp_freq_domain);
    }
}

template <unsigned long M>
void process_comp_unit(benchmark::State& state, float (&time_domain)[M]) {
    static_assert((M % 2) == 0, "Interleaved complex buffers must have even length.");
    constexpr unsigned long N = M / 2;
    alignas(MY_MAX_ALIGNMENT) float temp_time_domain[M];
    alignas(MY_MAX_ALIGNMENT) float temp_freq_domain[M] = {0};
    std::copy(time_domain, time_domain + M, temp_time_domain);

    for (auto _ : state) {
        if constexpr (N == 2) {
            FFT::fft_2(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 3) {
            FFT::fft_3(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 4) {
            FFT::fft_4(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 5) {
            FFT::fft_5(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 6) {
            FFT::fft_6(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 7) {
            FFT::fft_7(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 8) {
            FFT::fft_8(temp_time_domain, temp_freq_domain);
        }
    }
}

BENCHMARK_DEFINE_F(Case2, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case3, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case4, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case5, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case6, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case7, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }
BENCHMARK_DEFINE_F(Case8, CompUnitRaw)(benchmark::State& state) { process_comp_unit(state, time_domain); }

BENCHMARK_DEFINE_F(Case2, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case3, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case4, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case5, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case6, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case7, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(Case8, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }

BENCHMARK_DEFINE_F(PowerOf2, RawFFT)(benchmark::State& state) { process_fft(state, time_domain); }
BENCHMARK_DEFINE_F(PowerOf2, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }
BENCHMARK_DEFINE_F(PowerOf2, FFTWWithAlloc)(benchmark::State& state) { process_fftw_with_alloc(state, time_domain); }

BENCHMARK_DEFINE_F(MediumPrime, RawFFT)(benchmark::State& state) { process_fft(state, time_domain); }
BENCHMARK_DEFINE_F(MediumPrime, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }

BENCHMARK_DEFINE_F(MersennePrime, RawFFT)(benchmark::State& state) { process_fft(state, time_domain); }
BENCHMARK_DEFINE_F(MersennePrime, FFTW)(benchmark::State& state) { process_fftw(state, time_domain); }

BENCHMARK_REGISTER_F(Case2, CompUnitRaw);
BENCHMARK_REGISTER_F(Case2, FFTW);
BENCHMARK_REGISTER_F(Case3, CompUnitRaw);
BENCHMARK_REGISTER_F(Case3, FFTW);
BENCHMARK_REGISTER_F(Case4, CompUnitRaw);
BENCHMARK_REGISTER_F(Case4, FFTW);
BENCHMARK_REGISTER_F(Case5, CompUnitRaw);
BENCHMARK_REGISTER_F(Case5, FFTW);
BENCHMARK_REGISTER_F(Case6, CompUnitRaw);
BENCHMARK_REGISTER_F(Case6, FFTW);
BENCHMARK_REGISTER_F(Case7, CompUnitRaw);
BENCHMARK_REGISTER_F(Case7, FFTW);
BENCHMARK_REGISTER_F(Case8, CompUnitRaw);
BENCHMARK_REGISTER_F(Case8, FFTW);

BENCHMARK_REGISTER_F(PowerOf2, RawFFT);
BENCHMARK_REGISTER_F(PowerOf2, FFTW);
BENCHMARK_REGISTER_F(PowerOf2, FFTWWithAlloc);

BENCHMARK_REGISTER_F(MediumPrime, RawFFT);
BENCHMARK_REGISTER_F(MediumPrime, FFTW);

BENCHMARK_REGISTER_F(MersennePrime, RawFFT);
BENCHMARK_REGISTER_F(MersennePrime, FFTW);

BENCHMARK_MAIN();
