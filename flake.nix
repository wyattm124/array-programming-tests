{
  description = "C++ development environment with Clang";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs { inherit system; };
      # Use llvmPackages.stdenv for non-macOS platforms, default otherwise
      isMac = pkgs.lib.strings.hasSuffix "darwin" system;
      stdenv = if isMac
               then pkgs.stdenv
               else pkgs.llvmPackages.stdenv;
      # Determine architecture macro based on system
      archMacro = if pkgs.lib.strings.hasPrefix "aarch64" system || pkgs.lib.strings.hasPrefix "arm" system
                  then "ARM"
                  else if pkgs.lib.strings.hasPrefix "x86_64" system || pkgs.lib.strings.hasPrefix "i686" system
                  then "x86_64"
                  else "DEFAULT";
      doctest = stdenv.mkDerivation {
        pname = "doctest";
        version = "2.4.12";
        src = pkgs.fetchgit {
          url = "https://github.com/doctest/doctest.git";
          rev = "v2.4.12";
          sha256 = "sha256-Fxs1EWydhqN9whx+Cn4fnZ4fhCEQvFgL5e9TUiXlnq8="; 
        };
        buildInputs = [ pkgs.cmake pkgs.ninja ];
        CXXFLAGS = "-Wno-unsafe-buffer-usage";
        CMAKE_CXX_FLAGS = [ "-DDOCTEST_WITH_TESTS=OFF" "-DTREAT_WARNINGS_AS_ERRORS=OFF" ];
      };
      googlebench = stdenv.mkDerivation {
        pname = "google-benchmark";
        version = "1.9.4";
        src = pkgs.fetchgit {
          url = "https://github.com/google/benchmark.git";
          rev = "v1.9.4";
          sha256 = "sha256-P7wJcKkIBoWtN9FCRticpBzYbEZPq71a0iW/2oDTZRU=";
        };
        buildInputs = [ pkgs.cmake pkgs.ninja ];
        cmakeFlags = [
          "-DBENCHMARK_ENABLE_TESTING=OFF"
        ];
      };
    in rec {
      devShells.default = pkgs.mkShell {
        stdenv = stdenv;
        buildInputs = [
          stdenv.cc
          doctest
          googlebench
          pkgs.fftwFloat
        ] ++ pkgs.lib.optionals (!isMac) [
          pkgs.linuxPackages.perf
          pkgs.llvmPackages.llvm
        ];

        shellHook = ''
          unset NIX_ENFORCE_NO_NATIVE
          export CC=clang
          export CXX=clang++
          export CXXFLAGS="-DARCH_${archMacro} $CXXFLAGS"
          build_fft () {
            clang++ -std=c++23 -O3 -march=native test/fft_tests.cpp -lfftw3f -o bin/fft_tests && \
            clang++ -std=c++23 -O3 -march=native bench/fft_bench.cpp bench/fft_comp_unit.cpp -lbenchmark -pthread -lfftw3f -o bin/fft_bench && \
            clang++ -std=c++23 -O3 -g -march=native bench/fft_profile.cpp bench/fft_comp_unit.cpp -o bin/fft_profile
          }
          mca_timeline () {
            clang++ -std=c++23 -O3 -march=native bench/fft_profile.cpp -S -o - | llvm-mca -skip-unsupported-instructions=lack-sched --timeline
          } 
          build_tools () {
            clang++ -std=c++23 -O2 test/cheb_tests.cpp -o bin/cheb_tests && \
            clang++ -std=c++23 -O2 tools/smooth_dist.cpp -o bin/smooth_dist && \
            clang++ -std=c++23 -O2 tools/roots_printer.cpp -o bin/roots_printer
          }
        '' + pkgs.lib.optionalString (!isMac) ''
          perf_layer_analysis () {
            perf record -g --call-graph dwarf -e branch-misses,cache-misses,cycles ./bin/fft_profile && \
            perf report --stdio --symbol-filter=FFT::.*
          }
          perf_annotate_layer () {
            perf record -g --call-graph dwarf ./bin/fft_profile && \
            perf annotate --stdio
          }
        '';
      };

      packages.default = devShells.default; 
    });
}