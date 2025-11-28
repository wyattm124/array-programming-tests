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

Once Nix is installed, first build the development shell:

```bash
nix build
```

Then enter the shell with

```bash
nix develop
```

This will provide you with all the necessary dependencies including:
- doctest (testing framework)
- Google Benchmark
- FFTW (Fast Fourier Transform library)
- gperftools (performance profiling tools)
- Graphviz and gv (graph visualization tools)

The shell should not have to be rebuilt unless the nix flake is modified,
but the shell will have to be active to run any of the command aliases or
run any of the binaries.

## Shell Aliases

The development shell provides several convenient functions for
building and working with the project. Make sure to run these commands
from the root of the repository:

### `build_fft`

Builds all FFT-related executables:

- `bin/fft_tests` - FFT test suite
- `bin/fft_bench` - FFT benchmarking tool
- `bin/fft_profile` - FFT profiling tool

**Usage:**
```bash
build_fft
```

**Note:** These are compiled with `-O3` optimization level.

### `mca_timeline`

Runs LLVM Machine Code Analyzer (llvm-mca) with timeline analysis on
the FFT profiling code. This generates a detailed timeline view of
instruction scheduling and execution.

**Usage:**
```bash
mca_timeline
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

### `build_tools`

Builds utility tools:

- `bin/cheb_tests` - Chebyshev approximation tests
- `bin/smooth_dist` - Smooth distribution tool
- `bin/roots_printer` - Roots printing utility

**Usage:**
```bash
build_tools
```

**Note:** These are compiled with `-O2` optimization level.

These are assorted tools for exploring things such as the Chebyshev
polynomial approximation of functions, the number theoretic smoothness
of particular numbers as well as their distribution within specific
intervals as well as printing specific roots of unity with specified
precision.

## Project Structure

- `src/` - Source code headers (FFT, complex types, prime
  factorization)
- `test/` - Test files
- `bench/` - Benchmarking and profiling code
- `tools/` - Utility tools
- `bin/` - Compiled executables (created after building)

## Building Everything

To build all components:

```bash
build_fft
build_tools
```

