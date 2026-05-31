# Array Programming Tests

This repository contains C++ experiments for a form of array
programming that I have found helpful and worth exploring for common
numerical computing algorithms in areas such as computer vision and
signal processing. These experiments focus on doing a few practical
things well:

1. Aligning data arrays and types for easy translation to SIMD
   operations.
2. Preferring multiply and add over other mathematical operations
   wherever possible because these operations are easiest to reason about
   and the fastest for hardware to execute.
3. Compile time sizing of data arrays for predetermined memory
   requirements, predictable code the compiler can easily optimize,
   and effective compile time programming by an implementer. Many
   optimization techniques for this problem space rely on array sizes,
   so knowing array sizes at compile time becomes critical.
4. Easily interoperating with standard types, while offering optimized
   types that allow algorithms to operate at max speed.
5. Minimizing dependencies on the C++ standard library to simplify
   code, and allow complete independence if necessary.

And doing a few things algorithmically well:

1. Writing algorithms in array operation "layers" as you may see in an
   ONNX model because this representation is an effective and
   generalizable approach for these types of algorithms.
   This model also translates well to similar problem spaces like
   machine learning where ML models are commonly thought
   of in convolutional layers, or other similar "layers".
2. Writing transpositions as affine transformations wherever practical,
   as these forms are easiest to reason about with the polyhedral
   compilation model. Transpositions typically form most of the core
   operations in this problem space, so representing them in a
   understandable, composable, and easy to manipulate form is critical.
3. Specifying and computing as many algorithm details as possible at
   compile time. Compile time compute through constexpr and templates
   allows many aspects of an algorithm to be checked for correctness
   before runtime, and for runtime to require as little work as
   possible for peak performance. It also has the added benefit
   of compiling only the minimum code needed for a particular application.

This repository starts with implementing an FFT for a few reasons.
1. FFTs can be computed with well understood algorithms.
2. FFTs have many readily available and highly performant
   implementations, like FFTW3
3. FFTs are easy to define test cases and benchmarks for as their
   input is simply a summation of pure sine waves.
4. FFTs' highly symmetric nature requires one to experiment heavily
   with transpositions which I would argue are the key to optimizing
   many other similar problems.

The goal of this repository is to find and refine some techniques that
someone who is familiar with the pure math and algorithmic nature of a
problem in this space can leverage to write normal, portable, and
modern C++ (C++23 and newer) that operates in the ballpark of 85% of
a CPU's theoretically practical limits. In this way, such an
implementer can be confident in taking an idea from a math textbook
to a portable C++ implementation with great performance and
correctness.

Originally I thought pure C++ could get much closer to the greatest
possible performance, but after further investigation I realized that
much better performance is possible through carefully crafting assembly
or intrinsics based on the original mathematics of the problem and using
C++ only for the scaffolding of the overall solution. So, the original
goal has shifted some to include investigating what is required for, and
possible with, carefully crafted, and CPU specifc, code beyond
standard C++.

Note that this repo focuses on CPU computing, but the computing model
it focuses on is pretty general and could be used as a start for GPU
programming in this problem space as well.

## Prerequisites

### Installing Nix

This project uses [Nix](https://nixos.org/) for dependency management
and reproducible builds. To install Nix, you can use the [Determinate
Systems Nix Installer](https://github.com/DeterminateSystems/nix-installer),
which provides a fast and reliable way to install Nix with flakes
support.

## Getting Started

Once Nix is installed, enter the development shell from the repository root:

```bash
nix develop
```

The dev shell provides the compiler, Ninja, FFTW, Google Benchmark,
YAML tooling for generated sources, and a `pre-commit` hook that runs
`clang-format` on C and C++ files. It also exports architecture-specific
compile flags through `FFT_ARCH_FLAGS` so x86 builds use `AVX2`/`FMA`
and ARM builds use `NEON` when available.

From inside the shell, the normal workflow is Ninja-based:

```bash
ninja -f build/build.ninja build_fft
```

This builds the FFT executables and runs the FFT test suite as part of
the `build_fft` target.

The shell includes the main dependencies used by this repository:
- doctest (testing framework)
- Google Benchmark
- FFTW (Fast Fourier Transform library)
- yaml-cpp
- Python with PyYAML
- Ninja
- LLVM tools and Linux `perf` on supported Linux hosts

The shell does not define build aliases. Use the `build/build.ninja`
targets documented below.

## Build Commands

Run all commands below from the repository root.

### Build FFT targets

Build the FFT executables and run the FFT tests:

```bash
ninja -f build/build.ninja build_fft
```

This target builds:
- `bin/fft_tests`
- `bin/fft_bench`
- `bin/fft_profile`

### Run only the benchmark executable

```bash
ninja -f build/build.ninja run_fft_bench
```

To update benchmark baselines in `bench/specs/fft_bench.yaml` with the
measured CPU times from the current run:

```bash
FFT_BENCH_UPDATE_BASELINES=1 ninja -f build/build.ninja run_fft_bench
```

### Build utility targets

Build the non-FFT utilities and run `cheb_tests`:

```bash
ninja -f build/build.ninja build_tools
```

This target builds:
- `bin/cheb_tests`
- `bin/smooth_dist`
- `bin/roots_printer`

### Run all test suites

```bash
ninja -f build/build.ninja tests
```

This runs the FFT test suite and the utility test suite.

### Regenerate generated test and benchmark sources

The test and benchmark registration `.inc` files are generated from YAML
specs. To regenerate them explicitly:

```bash
ninja -f build/build.ninja generate_fft_tests
ninja -f build/build.ninja generate_fft_benchmarks
```

### `mca_timeline`

Runs LLVM Machine Code Analyzer (llvm-mca) with timeline analysis on
the FFT profiling code. This generates a detailed timeline view of
instruction scheduling and execution.

```bash
ninja -f build/build.ninja mca_timeline
```

This compiles `bench/fft_profile.cpp` to assembly and analyzes it
with `llvm-mca`, providing insights into CPU pipeline behavior and
instruction-level parallelism.

To select a section of code to observe with the tool, first modify
`bench/fft_profile.cpp` as necessary to make sure the code section
of interest will be compiled. Then put the `MCA_START` macro at the
beginning of the code selection, and `MCA_END` at the ending of the
code selection. Templated code sections may not work with these macros
as they cannot be nested, and are best if only defined once. These
macros are defined at the top of the `src/fft.hpp` file.

### `perf_layer_analysis`

On supported Linux hosts, collect a `perf report` focused on FFT layer
symbols:

```bash
ninja -f build/build.ninja perf_layer_analysis
```

### `perf_annotate_layer`

On supported Linux hosts, collect annotated `perf` output:

```bash
ninja -f build/build.ninja perf_annotate_layer
```

### Build everything

```bash
ninja -f build/build.ninja all
```

## Project Structure

- `src/` - FFT implementation headers, architecture detection, prime
  factorization helpers, and DFT layer specializations
- `src/fft.hpp` - main FFT plan and recursive decomposition logic
- `src/dft_default.h` - default DFT layer implementations used when no
  SIMD specialization exists
- `src/dft_AVX2.h` - AVX2-specific DFT layer specializations
- `test/` - doctest-based test sources
- `test/specs/` - YAML specifications used to generate FFT test cases
- `test/generated/` - generated `.inc` files included by `test/fft_tests.cpp`
- `bench/` - benchmark and profiling sources
- `bench/specs/fft_bench.yaml` - YAML benchmark specification and stored
  baseline CPU times
- `bench/generated/` - generated benchmark registration `.inc` files
- `tools/` - helper programs and code generators
- `tools/generate_fft_doctest.py` - generates FFT doctest includes from
  `test/specs/*.yaml`
- `tools/generate_fft_benchmarks.py` - generates benchmark registration
  includes from `bench/specs/fft_bench.yaml`
- `build/build.ninja` - canonical local build graph used by the repo
- `bin/` - linked executables produced by Ninja
- `build/obj/` - intermediate object files produced by Ninja

## Spec-Driven Generation

FFT tests and FFT benchmark registrations are not written entirely by
hand. Instead, the repo uses YAML specs as the source of truth:

- `test/specs/fft.yaml` controls generated FFT doctest cases and their
  tolerances for each SIMD mode
- `bench/specs/fft_bench.yaml` controls generated benchmark registrations
  and stores benchmark baselines used for delta reporting

When these spec files change, Ninja regenerates the matching files in
`test/generated/` and `bench/generated/` before compiling the owning
targets.

## Common Outputs

- `bin/fft_tests` - doctest-based FFT validation suite
- `bin/fft_bench` - Google Benchmark runner with spec-backed baseline reporting
- `bin/fft_profile` - focused profiling entry point for LLVM MCA and `perf`
- `bin/cheb_tests` - tests for Chebyshev approximation helpers
- `bin/smooth_dist` - smooth-number exploration tool
- `bin/roots_printer` - roots-of-unity printing utility

## Building Everything

```bash
ninja -f build/build.ninja all
```
