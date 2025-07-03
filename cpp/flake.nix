{
  description = "C++ development environment with Clang";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }: {
    devShell.defaultPackage = nixpkgs.mkShell {
      buildInputs = [
        nixpkgs.clang
        nixpkgs.cmake
        nixpkgs.gnumake
        nixpkgs.git
        nixpkgs.ninja
      ];

      shellHook = ''
        export CC=clang
        export CXX=clang++
      ''; 
    };

    packages.default = nixpkgs.mkShell {
      buildInputs = [
        nixpkgs.clang
        nixpkgs.cmake
      ];
    };
  };
}