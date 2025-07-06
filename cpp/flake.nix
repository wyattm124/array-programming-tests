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
        buildInputs = [ pkgs.cmake pkgs.ninja pkgs.libcxx ];
        CXXFLAGS = "-Wno-unsafe-buffer-usage";
        CMAKE_CXX_FLAGS = [ "-DDOCTEST_WITH_TESTS=OFF" "-DTREAT_WARNINGS_AS_ERRORS=OFF" ];
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
          pkgs.libcxx
        ];

        shellHook = ''
          export CC=clang
          export CXX=clang++
        '';
      };

      packages.default = pkgs.mkShell {
        buildInputs = [
          pkgs.clang
          pkgs.cmake
          doctest
          pkgs.libcxx
        ];
      };
    });
}