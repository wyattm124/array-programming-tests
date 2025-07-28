{
  description = "C++ development environment with Clang";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs { inherit system; };
      doctest = pkgs.stdenv.mkDerivation {
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
      googlebench = pkgs.stdenv.mkDerivation {
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
    in {
      devShells.default = pkgs.mkShell {
        buildInputs = [
          pkgs.clang
          pkgs.cmake
          pkgs.gnumake
          pkgs.git
          pkgs.ninja
          doctest
          googlebench
          pkgs.libcxx
          pkgs.fftwFloat
          pkgs.gperftools
          pkgs.graphviz
          pkgs.gv
          pkgs.llvmPackages_latest.llvm
        ];

        shellHook = ''
          export CC=clang
          export CXX=clang++
          build () {
            clang++ -std=c++23 -O3 fft_tests.cpp -o ../bin/fft_tests && \
            clang++ -std=c++23 -O3 fft_bench.cpp -lbenchmark -pthread -lfftw3f -o ../bin/fft_bench && \
            clang++ -std=c++23 -O3 fft_profile.cpp -lprofiler -o ../bin/fft_profile && \
            clang++ -std=c++23 -O3 cheb_tests.cpp -o ../bin/cheb_tests && \
            clang++ -std=c++23 -O3 prime_factor.cpp -o ../bin/prime_factor
          }
          mca_timeline () {
            clang++ fft_profile.cpp -O2 -S -o - | llvm-mca -skip-unsupported-instructions=lack-sched --timeline
          }
        '';
      };

      packages.default = pkgs.mkShell {
        buildInputs = [
          pkgs.clang
          pkgs.cmake
          doctest
          googlebench
          pkgs.libcxx
          pkgs.fftwFloat
          pkgs.gperftools
          pkgs.graphviz
          pkgs.gv
          pkgs.llvmPackages_latest.llvm
        ];
      };
    });
}