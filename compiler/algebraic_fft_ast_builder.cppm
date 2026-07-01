export module fft.compiler.graph;

import <memory>;
import <string>;
import <utility>;
import <vector>;

import fft.compiler.prime_factor;

export namespace fft_compiler {

struct Value {
  virtual ~Value() = default;
  [[nodiscard]] virtual std::string DotLabel() const = 0;
};

struct Zero : Value {
  [[nodiscard]] std::string DotLabel() const override { return "0"; }
};

struct RootOfUnity : Value {
  RootOfUnity(unsigned int root_num, unsigned int root_denom)
      : num(root_num), denom(root_denom) {}

  [[nodiscard]] std::string DotLabel() const override {
    return "ω^(" + std::to_string(num) + "/" + std::to_string(denom) + ")";
  }

  unsigned int num = 0;
  unsigned int denom = 1;
};

struct Input : Value {
  explicit Input(unsigned int input_index) : index(input_index) {}

  [[nodiscard]] std::string DotLabel() const override {
    return "x_" + std::to_string(index);
  }

  unsigned int index = 0;
};

struct AlgebraicASTNode {
  virtual ~AlgebraicASTNode() = default;
  [[nodiscard]] virtual std::string DotLabel() const = 0;
  [[nodiscard]] virtual std::string DotShape() const { return "ellipse"; }
  [[nodiscard]] virtual std::vector<std::shared_ptr<AlgebraicASTNode>>
  DotChildren() const {
    return {};
  }
};

struct ValueASTNode : AlgebraicASTNode {
  ValueASTNode(std::shared_ptr<Value> first_value,
               std::shared_ptr<Value> last_value)
      : first(std::move(first_value)), last(std::move(last_value)) {}

  [[nodiscard]] std::string DotLabel() const override {
    return "(" + first->DotLabel() + ", " + last->DotLabel() + ")";
  }

  [[nodiscard]] std::string DotShape() const override { return "box"; }

  std::shared_ptr<Value> first;
  std::shared_ptr<Value> last;
};

struct BinaryASTNode : AlgebraicASTNode {
  BinaryASTNode(std::shared_ptr<AlgebraicASTNode> left,
                std::shared_ptr<AlgebraicASTNode> right)
      : lhs(std::move(left)), rhs(std::move(right)) {}

  [[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  DotChildren() const override {
    return {lhs, rhs};
  }

  std::shared_ptr<AlgebraicASTNode> lhs;
  std::shared_ptr<AlgebraicASTNode> rhs;
};

struct Mult : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "*"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct Add : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "+"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct PluckOuter : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "outer"; }
};

struct PluckInner : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "inner"; }
};

struct PluckFirst : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "first"; }
};

struct PluckLast : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "last"; }
};

struct AlgebraicDFTKernelGraph {
  unsigned int dft_size = 0;
  unsigned int dfts_per_kernel = 1;
  std::vector<std::pair<unsigned int, unsigned int>> prime_factor_exponents;
  std::vector<std::shared_ptr<AlgebraicASTNode>> nodes;
};

class AlgebraicDFTKernelGraphBuilder {
public:
  [[nodiscard]] AlgebraicDFTKernelGraph Build(
      unsigned int dft_size, unsigned int dfts_per_kernel) const {
    AlgebraicDFTKernelGraph graph;
    graph.dft_size = dft_size;
    graph.dfts_per_kernel = dfts_per_kernel;
    graph.prime_factor_exponents =
        prime_factor::get_prime_factor_powers(dft_size);
    graph.nodes = genGraph(dft_size, dfts_per_kernel);
    return graph;
  }

private:
  [[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  genGraph(unsigned int dft_size, unsigned int dfts_per_kernel) const {
    std::vector<std::shared_ptr<AlgebraicASTNode>> roots;
    if (dft_size == 0) {
      return roots;
    }

    const auto input_nodes = genInputNodes(dft_size);
    const auto shuffled_inputs = shuffleInputNodes(input_nodes, dft_size);
    const auto tree_sum = treeReduceSum(shuffled_inputs, 1, 0);

    for (unsigned int kernel_index = 0; kernel_index < dfts_per_kernel;
         ++kernel_index) {
      if (tree_sum) {
        roots.push_back(tree_sum);
      }

      if (kernel_index + 1 < dfts_per_kernel && dft_size > 1) {
        auto kernel_root = std::make_shared<RootOfUnity>(kernel_index + 1, dft_size);
        auto kernel_tag = std::make_shared<ValueASTNode>(kernel_root, std::make_shared<Zero>());
        roots.push_back(kernel_tag);
      }
    }

    return roots;
  }

  [[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  genInputNodes(unsigned int dft_size) const {
    std::vector<std::shared_ptr<AlgebraicASTNode>> nodes;
    for (unsigned int i = 0; i + 1 < dft_size; i += 2) {
      nodes.emplace_back(std::make_shared<ValueASTNode>(
          std::make_shared<Input>(i), std::make_shared<Input>(i + 1)));
    }

    if (dft_size % 2 != 0) {
      nodes.emplace_back(std::make_shared<ValueASTNode>(
          std::make_shared<Input>(dft_size - 1), std::make_shared<Zero>()));
    }

    return nodes;
  }

  [[nodiscard]] static std::vector<std::shared_ptr<AlgebraicASTNode>>
  shuffleInputNodes(const std::vector<std::shared_ptr<AlgebraicASTNode>> &input_layer,
                    unsigned int dft_size) {
    std::vector<std::shared_ptr<AlgebraicASTNode>> result;
    const unsigned int len = static_cast<unsigned int>(input_layer.size());
    if (len == 0) {
      return result;
    }

    unsigned int left = 0;
    unsigned int right = len - 1;

    if (dft_size % 2) {
      while (left < right) {
        result.emplace_back(
            std::make_shared<PluckOuter>(input_layer[left], input_layer[right]));
        result.emplace_back(
            std::make_shared<PluckInner>(input_layer[left], input_layer[right]));
        ++left;
        --right;
      }
      if (left == right) {
        result.push_back(input_layer[left]);
      }
      return result;
    }

    if (len % 2) {
      result.emplace_back(std::make_shared<PluckOuter>(input_layer[0],
                                                       input_layer[len / 2]));
    } else {
      result.emplace_back(std::make_shared<PluckFirst>(input_layer[0],
                                                       input_layer[len / 2]));
    }

    bool take_right = true;
    while (left < right) {
      if (take_right) {
        result.emplace_back(std::make_shared<PluckLast>(input_layer[left++],
                                                        input_layer[right]));
      } else {
        result.emplace_back(std::make_shared<PluckFirst>(input_layer[left],
                                                         input_layer[right--]));
      }
      take_right = !take_right;
    }

    return result;
  }

  [[nodiscard]] static std::shared_ptr<AlgebraicASTNode>
  treeReduceSum(const std::vector<std::shared_ptr<AlgebraicASTNode>> &input,
                unsigned int stride, unsigned int offset) {
    if (!stride) {
      return {};
    }

    std::vector<std::shared_ptr<AlgebraicASTNode>> current_layer;
    for (unsigned int i = offset; i < input.size(); i += stride) {
      current_layer.push_back(input[i]);
    }

    if (current_layer.empty()) {
      return {};
    }

    while (current_layer.size() > 1) {
      std::vector<std::shared_ptr<AlgebraicASTNode>> next_layer;
      for (std::size_t i = 0; i + 1 < current_layer.size(); i += 2) {
        next_layer.emplace_back(
            std::make_shared<Add>(current_layer[i], current_layer[i + 1]));
      }

      if (current_layer.size() % 2) {
        next_layer.push_back(current_layer.back());
      }

      current_layer = std::move(next_layer);
    }

    return current_layer.front();
  }
};

using FFTKernelGraphBuilder = AlgebraicDFTKernelGraphBuilder;
using InstructionGraph = AlgebraicDFTKernelGraph;

} // namespace fft_compiler
