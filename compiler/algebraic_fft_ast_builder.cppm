export module fft.compiler.graph;

import <memory>;
import <utility>;
import <vector>;

import fft.compiler.prime_factor;

export namespace fft_compiler {

struct AlgebraicASTNode {
  virtual ~AlgebraicASTNode() = default;
};

struct InstructionGraph {
  unsigned int dft_size = 0;
  unsigned int dfts_per_kernel = 1;
  std::vector<std::pair<unsigned int, unsigned int>> prime_factor_exponents;
  std::vector<std::shared_ptr<AlgebraicASTNode>> nodes;
};

class FFTKernelGraphBuilder {
public:
  [[nodiscard]] InstructionGraph Build(unsigned int dft_size,
                                       unsigned int dfts_per_kernel) const {
    InstructionGraph graph;
    graph.dft_size = dft_size;
    graph.dfts_per_kernel = dfts_per_kernel;
    graph.prime_factor_exponents =
        prime_factor::get_prime_factor_powers(dft_size);

    // TODO: Build the algebraic FFT AST and populate graph.nodes.

    return graph;
  }
};

using AlgebraicDFTKernelGraphBuilder = FFTKernelGraphBuilder;

} // namespace fft_compiler
