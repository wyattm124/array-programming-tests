{
  description = "C++ development environment with Clang";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    git-hooks.url = "github:cachix/git-hooks.nix";
  };

  outputs = { self, nixpkgs, flake-utils, git-hooks }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs { inherit system; };
      # Use llvmPackages.stdenv for non-macOS platforms, default otherwise
      isMac = pkgs.lib.strings.hasSuffix "darwin" system;
      stdenv = if isMac
               then pkgs.stdenv
               else pkgs.llvmPackages.stdenv;
      archCompileFlags =
        if pkgs.lib.strings.hasPrefix "aarch64" system || pkgs.lib.strings.hasPrefix "arm64" system
        then "-march=armv8-a+simd"
        else if pkgs.lib.strings.hasPrefix "arm" system
        then "-mfpu=neon"
        else if pkgs.lib.strings.hasPrefix "x86_64" system || pkgs.lib.strings.hasPrefix "i686" system
        then "-mavx2 -mfma"
        else "";
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
      pythonWithYaml = pkgs.python3.withPackages (ps: [ ps.pyyaml ]);
      pre-commit-check = git-hooks.lib.${system}.run {
        src = ./.;
        hooks = {
          clang-format = {
            enable = true;
            files = "\\.(cpp|cc|cxx|c|hpp|hh|hxx|h)$";
          };
        };
      };
    in rec {
      devShells.default = pkgs.mkShell {
        stdenv = stdenv;
        buildInputs = [
          pkgs.ninja
          stdenv.cc
          doctest
          googlebench
          pkgs.fftwFloat
          pkgs.yaml-cpp
          pythonWithYaml
        ] ++ pre-commit-check.enabledPackages
          ++ pkgs.lib.optionals (!isMac) [
          pkgs.linuxPackages.perf
          pkgs.llvmPackages.llvm
        ];

        shellHook = ''
          ${pre-commit-check.shellHook}
          unset NIX_ENFORCE_NO_NATIVE
          export CC=clang
          export CXX=clang++
          export FFT_ARCH_FLAGS="${archCompileFlags}"
        '';
      };

      packages.default = devShells.default; 
    });
}
