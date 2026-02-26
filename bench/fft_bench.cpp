#include "fft_comp_unit.hpp"
#include "../src/fft.hpp"
#include <benchmark/benchmark.h>
#include <fftw3.h>

#include <vector>

/// Base Class for Sizes
template <unsigned int M>
class Case : public benchmark::Fixture {
public:
    constexpr static unsigned int N = M;
protected:
    alignas(MY_MAX_ALIGNMENT) std::complex<float> time_domain[N];
    alignas(MY_MAX_ALIGNMENT) std::complex<float> freq_domain[N];
    
    void SetUp(const benchmark::State& state) {
        for (unsigned int i = 0; i < N; i++) {
            time_domain[i] = {0, 0};
            freq_domain[i] = {0, 0};
        }
        FFT::wave_gen_lcg(time_domain, freq_domain, N);
    }
};

// Base Cases
using Case2 = Case<2>;
using Case3 = Case<3>;
using Case4 = Case<4>;
using Case5 = Case<5>;
using Case6 = Case<6>;
using Case7 = Case<7>;
using Case8 = Case<8>;

// Power of 2 cases
using PowerOf2 = Case<8192>;

// Medium Prime cases
using MediumPrime = Case<53>;

// Mersenne Prime cases
using MersennePrime = Case<8191>;


/// Base Classes for Processing Flows
template <unsigned long N>
void process_fftw(benchmark::State& state, std::complex<float> (&time_domain)[N]) {
    // Configure input with waves
    fftwf_complex *out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
    fftwf_plan p;

    // Set the plan
    p = fftwf_plan_dft_1d(N,
        reinterpret_cast<fftwf_complex*>(time_domain),
        out, FFTW_FORWARD, FFTW_MEASURE);

    // Benchmark
    for (auto _ : state) {
        fftwf_execute(p);
    }

    // Cleanup
    fftwf_destroy_plan(p);
    fftwf_free(out);
}

template <unsigned long N>
void process_fftw_with_alloc(benchmark::State& state, std::complex<float> (&time_domain)[N]) {
    // Configure input with waves
    fftwf_complex *in, *out;
    fftwf_plan p;
    std::vector<fftwf_complex *> buffs;

    // Benchmark
    for (auto _ : state) {
        // Allocate fftw specific memory
        in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
        out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
        
        // Set the plan
        p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        // Copy the input to the appropriate location
        std::copy(
            reinterpret_cast<const float*>(time_domain),
            reinterpret_cast<const float*>(time_domain + N),
            reinterpret_cast<float*>(in)
        );
        
        // Execute the fft itself
        fftwf_execute(p);

        // Free memory later
        buffs.push_back(in);
        buffs.push_back(out);

        // Cleanup plan for reset
        fftwf_destroy_plan(p);
    }

    for (auto buf : buffs) {
        fftwf_free(buf);
    }
}

#ifdef ARCH_ARM
template <unsigned long N>
void process_neon(benchmark::State& state, std::array<std::complex<float>, N>& time_domain) {
    // Need to Copy Data to Correct Type
    std::array<float32x2_t, N> temp_time_domain;
    std::array<float32x2_t, N> temp_freq_domain;
    for (unsigned int i = 0; i < N; i++) {
        temp_time_domain[i] = {time_domain[i].real(), time_domain[i].imag()};
    }

    // Make sure to do initialization
    FFT::FFTPlan<float32x2_t>::Init<N>();

    for (auto _ : state) {
        FFT::FFTPlan<float32x2_t>::fft<N>(temp_time_domain.data(), temp_freq_domain.data());
    }
}
#endif

template <unsigned long N>
void process_stdcomplexwrap(benchmark::State& state, std::complex<float> (&time_domain)[N]) {
    FFT::StdComplexWrap<float> temp_freq_domain[N];

    // Make sure to do initialization
    FFT::FFTPlan<FFT::StdComplexWrap<float>>::Init<N>();

    for (auto _ : state) {
        FFT::FFTPlan<FFT::StdComplexWrap<float>>::fft<N>(
            static_cast<FFT::StdComplexWrap<float>*>(time_domain),
            temp_freq_domain
        );
    }
}

template <unsigned long N>
void process_customcomplex(benchmark::State& state, std::complex<float> (&time_domain)[N]) {
    // Need to Copy Data to Correct Type
    alignas(MY_MAX_ALIGNMENT) FFT::Complex temp_time_domain[N];
    alignas(MY_MAX_ALIGNMENT) FFT::Complex temp_freq_domain[N];
    for (unsigned int i = 0; i < N; i++) {
        temp_time_domain[i] = {time_domain[i].real(), time_domain[i].imag()};
    }

    // Make sure to do initialization
    FFT::FFTPlan<FFT::Complex>::Init<N>();

    for (auto _ : state) {
        FFT::FFTPlan<FFT::Complex>::fft<N>(temp_time_domain, temp_freq_domain);
    }
}

template <unsigned long N>
void process_comp_unit(benchmark::State& state, std::complex<float> (&time_domain)[N]) {
    // Need to Copy Data to Correct Type
    alignas(MY_MAX_ALIGNMENT) FFT::Complex temp_time_domain[N];
    alignas(MY_MAX_ALIGNMENT) FFT::Complex temp_freq_domain[N];
    for (unsigned int i = 0; i < N; i++) {
        temp_time_domain[i] = {time_domain[i].real(), time_domain[i].imag()};
    }

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

BENCHMARK_DEFINE_F(Case2, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case3, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case4, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case5, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case6, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case7, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case8, CompUnitComplex)(benchmark::State& state) {
    process_comp_unit(state, time_domain);
}

BENCHMARK_DEFINE_F(Case2, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case3, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case4, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case5, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case6, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case7, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(Case8, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

#ifdef ARCH_ARM
BENCHMARK_DEFINE_F(PowerOf2, Neon)(benchmark::State& state) {
    process_neon(state, time_domain);
}
#endif

BENCHMARK_DEFINE_F(PowerOf2, CustomComplex)(benchmark::State& state) {
    process_customcomplex(state, time_domain);
}

BENCHMARK_DEFINE_F(PowerOf2, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(PowerOf2, FFTWWithAlloc)(benchmark::State& state) {
    process_fftw_with_alloc(state, time_domain);
}

BENCHMARK_DEFINE_F(PowerOf2, StdComplexWrap)(benchmark::State& state) {
    process_stdcomplexwrap(state, time_domain);
}

BENCHMARK_DEFINE_F(MediumPrime, CustomComplex)(benchmark::State& state) {
    process_customcomplex(state, time_domain);
}

BENCHMARK_DEFINE_F(MediumPrime, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

BENCHMARK_DEFINE_F(MersennePrime, CustomComplex)(benchmark::State& state) {
    process_customcomplex(state, time_domain);
}

BENCHMARK_DEFINE_F(MersennePrime, FFTW)(benchmark::State& state) {
    process_fftw(state, time_domain);
}

// Base Case Benchmarks
BENCHMARK_REGISTER_F(Case2, CompUnitComplex);
BENCHMARK_REGISTER_F(Case2, FFTW);
BENCHMARK_REGISTER_F(Case3, CompUnitComplex);
BENCHMARK_REGISTER_F(Case3, FFTW);
BENCHMARK_REGISTER_F(Case4, CompUnitComplex);
BENCHMARK_REGISTER_F(Case4, FFTW);
BENCHMARK_REGISTER_F(Case5, CompUnitComplex);
BENCHMARK_REGISTER_F(Case5, FFTW);
BENCHMARK_REGISTER_F(Case6, CompUnitComplex);
BENCHMARK_REGISTER_F(Case6, FFTW);
BENCHMARK_REGISTER_F(Case7, CompUnitComplex);
BENCHMARK_REGISTER_F(Case7, FFTW);
BENCHMARK_REGISTER_F(Case8, CompUnitComplex);
BENCHMARK_REGISTER_F(Case8, FFTW);

// Power of 2 Benchmarks
#ifdef ARCH_ARM
BENCHMARK_REGISTER_F(PowerOf2, Neon);
#endif
BENCHMARK_REGISTER_F(PowerOf2, CustomComplex);
BENCHMARK_REGISTER_F(PowerOf2, FFTW);
BENCHMARK_REGISTER_F(PowerOf2, FFTWWithAlloc);
BENCHMARK_REGISTER_F(PowerOf2, StdComplexWrap);

// Medium Prime Benchmarks
BENCHMARK_REGISTER_F(MediumPrime, CustomComplex);
BENCHMARK_REGISTER_F(MediumPrime, FFTW);

// Mersenne Prime Benchmarks
BENCHMARK_REGISTER_F(MersennePrime, CustomComplex);
BENCHMARK_REGISTER_F(MersennePrime, FFTW);

// Template out main
BENCHMARK_MAIN();