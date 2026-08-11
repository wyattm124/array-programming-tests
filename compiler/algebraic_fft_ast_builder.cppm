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

struct UnaryASTNode : AlgebraicASTNode {
  UnaryASTNode(std::shared_ptr<AlgebraicASTNode> input)
      : in(std::move(input)){}

  [[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  DotChildren() const override {
    return {in};
  }

  std::shared_ptr<AlgebraicASTNode> in;
};

struct Mult : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "*"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct MultConj : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "*~"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct Add : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "+"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct Sub : BinaryASTNode {
  using BinaryASTNode::BinaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "-"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct FoldAdd : UnaryASTNode {
  using UnaryASTNode::UnaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "+"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct FoldSub : UnaryASTNode {
  using UnaryASTNode::UnaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "-"; }
  [[nodiscard]] std::string DotShape() const override { return "circle"; }
};

struct Negate : UnaryASTNode {
  using UnaryASTNode::UnaryASTNode;
  [[nodiscard]] std::string DotLabel() const override { return "Neg"; }
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

struct JumpListEntry {
  bool visited = false;
  std::vector<unsigned int> next_jump_factors = {};
};

[[nodiscard]] std::shared_ptr<AlgebraicASTNode>
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

[[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  genInputNodes(unsigned int dft_size) {
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

[[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
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

std::vector<JumpListEntry> createJumpList(unsigned int len) {
  std::vector<JumpListEntry> jump_list(len);
  for (unsigned int i = 0; i < len; i++) {
    unsigned int num = len - i;
    
    // Find the preceding stride that diffs by the smallest factor.
    // This is equivalent to finding a preceding stride with the greatest
    // GCD.
    unsigned int prime_factor = prime_factor::get_prime_factor(num);
    if (prime_factor > 1 && (num % prime_factor) == 0) {
      jump_list[(num/prime_factor) - 1]
        .next_jump_factors.push_back(prime_factor);
    }
  }
  return jump_list;
}

void createSumsListHelper(
  const std::vector<std::shared_ptr<AlgebraicASTNode>> &source_list,
  std::vector<std::vector<std::shared_ptr<AlgebraicASTNode>>> &sums_list,
  std::vector<JumpListEntry> &jump_list,
  unsigned int curr_val,
  unsigned int dft_size) {

  // Populate the sums list for this index
  std::vector<std::shared_ptr<AlgebraicASTNode>> sink_list;
  if (dft_size % curr_val) {
    sink_list = source_list;
  } else {
    unsigned int bin_size = dft_size/curr_val;
    for (unsigned int j = 0; j < bin_size && j < source_list.size(); j++) {
      if (j + bin_size < source_list.size()) {
        sink_list.push_back(treeReduceSum(source_list, bin_size, j));
      } else {
        sink_list.push_back(source_list[j]);
      }
    }
  }
  sums_list[curr_val] = sink_list;
  
  // Start DFS of all sums lists that use the current one
  //  based off of a jump factor
  auto &jump_list_entry = jump_list[curr_val - 1];
  for (const auto &jump_factor : jump_list_entry.next_jump_factors) { 
    
    createSumsListHelper(sink_list, sums_list, jump_list,
      curr_val * jump_factor, dft_size);
  }

  // Mark this jump as visited
  jump_list_entry.visited = true;
}

std::vector<std::vector<std::shared_ptr<AlgebraicASTNode>>> createSumsList(
  const std::vector<std::shared_ptr<AlgebraicASTNode>> &shuffled_input,
  unsigned int dft_size) {
  
  // Create the jump list with outputs of shared sums.
  //  There should be an entry for every input pair not starting with the
  //  0th scalar input. Note the jump list will be 1 shorter than the sums
  //  lists.
  auto jump_list = createJumpList(shuffled_input.size() - 1);

  // The first output is always the sum of the inputs, and there will be an
  //  and sum list for every pair of inputs.
  std::vector<std::vector<std::shared_ptr<AlgebraicASTNode>>> sums_lists(shuffled_input.size());
  sums_lists[0] = {{shuffled_input}};

  // Go through the jump list to create the shared sums
  for (unsigned int i = 0; i < jump_list.size(); i++) {
    
    // Skip this entry if it has already been visited.
    //  This should only be the case for indexes relatively prime
    //  to the DFT size.
    if (jump_list[i].visited)
      continue; 

    // DFS through next jumps starting at this index
    createSumsListHelper(shuffled_input, sums_lists, jump_list, i + 1, dft_size);
  }

  return sums_lists;
}

std::vector<std::shared_ptr<AlgebraicASTNode>> sumsListDotProd(
  const std::vector<std::vector<std::shared_ptr<AlgebraicASTNode>>> &sums_list,
  unsigned int dft_size) {

  // Populate the results that are strictly sums
  std::vector<std::shared_ptr<AlgebraicASTNode>> result(dft_size);
  result[0] = treeReduceSum(sums_list[0], 1, 0);
  if ((dft_size - 1) % 2) {
    result[dft_size/2] = result[0];
  }

  // Keep the roots cached to avoid re-calculating them 
  std::vector<std::shared_ptr<AlgebraicASTNode>> roots;
  for (unsigned int i = 1; i < (dft_size + 1)/2; i++) {
    unsigned int conj_index = dft_size - i;
    roots.push_back(std::make_shared<ValueASTNode>(
      std::make_shared<RootOfUnity>(i, dft_size),
      std::make_shared<RootOfUnity>(conj_index, dft_size)
    ));
  }
 
  for (unsigned int i = 1; i < (dft_size - i); i++) {
    std::vector<std::shared_ptr<AlgebraicASTNode>> prod_map{sums_list[i][0]};
    std::vector<std::shared_ptr<AlgebraicASTNode>> prod_map_conj{sums_list[i][0]};

    for (unsigned int j = 1; j < sums_list[i].size(); j++) {
      unsigned int root_val = (i * j) % dft_size;

      if (root_val <= roots.size()) {
        prod_map.push_back(std::make_shared<Mult>(sums_list[i][j], roots[root_val - 1]));
        prod_map_conj.push_back(std::make_shared<MultConj>(sums_list[i][j], roots[root_val - 1]));
      } else if (root_val == roots.size() + 1) {
        auto neg_input = std::make_shared<Negate>(sums_list[i][j]);
        prod_map.push_back(neg_input);
        prod_map_conj.push_back(neg_input);
      } else {
        root_val = (2 * roots.size() + 1) - root_val;
        prod_map.push_back(std::make_shared<MultConj>(sums_list[i][j], roots[root_val]));
        prod_map_conj.push_back(std::make_shared<Mult>(sums_list[i][j], roots[root_val]));
      }
    }
    
    result[i] = treeReduceSum(prod_map, 1, 0);
    result[dft_size - i] = treeReduceSum(prod_map_conj, 1, 0);
  }
 
  return result;
}

std::vector<std::shared_ptr<AlgebraicASTNode>> foldToFinalResult(
  const std::vector<std::shared_ptr<AlgebraicASTNode>> &result_pairs_list, 
  unsigned int dft_size) {
  std::vector<std::shared_ptr<AlgebraicASTNode>> result(
    result_pairs_list.size());
  result[0] = std::make_shared<FoldAdd>(result_pairs_list[0]);
  for (unsigned int i = 1; i < (dft_size - i); i++) {
    result[i] = std::make_shared<FoldAdd>(result_pairs_list[i]);
    result[dft_size - i] = std::make_shared<FoldAdd>(result_pairs_list[dft_size - i]);
  }

  if ((dft_size - 1) % 2) {
    result[dft_size/2] = std::make_shared<FoldSub>(result_pairs_list[dft_size/2]);
  }
  return result;
}

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
    graph.nodes = genGraph(dft_size);
    return graph;
  }

private:
  [[nodiscard]] std::vector<std::shared_ptr<AlgebraicASTNode>>
  genGraph(unsigned int dft_size) const {
    std::vector<std::shared_ptr<AlgebraicASTNode>> roots;
    if (dft_size == 0) {
      return roots;
    }

    const auto input_nodes = genInputNodes(dft_size);
    const auto shuffled_inputs = shuffleInputNodes(input_nodes, dft_size);
    const auto sums_list = createSumsList(shuffled_inputs, dft_size);
    const auto paired_results = sumsListDotProd(sums_list, dft_size);
    const auto final_results = foldToFinalResult(paired_results, dft_size); 

    return final_results;
  } 
};

using FFTKernelGraphBuilder = AlgebraicDFTKernelGraphBuilder;
using InstructionGraph = AlgebraicDFTKernelGraph;

} // namespace fft_compiler
