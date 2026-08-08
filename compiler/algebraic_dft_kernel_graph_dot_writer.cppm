export module fft.compiler.dot;

import <filesystem>;
import <fstream>;
import <memory>;
import <sstream>;
import <string>;
import <string_view>;
import <vector>;

import fft.compiler.graph;

namespace {

std::string escape_dot_label(std::string_view text) {
  std::string escaped;
  escaped.reserve(text.size());

  for (const char ch : text) {
    switch (ch) {
    case '\\':
      escaped += "\\\\";
      break;
    case '"':
      escaped += "\\\"";
      break;
    case '\n':
      escaped += "\\n";
      break;
    default:
      escaped += ch;
      break;
    }
  }

  return escaped;
}

std::string ast_node_id(const std::shared_ptr<fft_compiler::AlgebraicASTNode> &node) {
  std::ostringstream id;
  id << "ast_" << static_cast<const void *>(node.get());
  return id.str();
}

bool was_visited(const std::vector<const fft_compiler::AlgebraicASTNode *> &visited,
                 const fft_compiler::AlgebraicASTNode *node) {
  for (const auto *seen : visited) {
    if (seen == node) {
      return true;
    }
  }
  return false;
}

void append_ast_subgraph(
    std::ostringstream &dot,
    const std::shared_ptr<fft_compiler::AlgebraicASTNode> &node,
    std::vector<const fft_compiler::AlgebraicASTNode *> &visited,
    std::vector<std::string> &leaf_node_ids) {
  if (!node || was_visited(visited, node.get())) {
    return;
  }

  visited.push_back(node.get());
  const std::string node_id = ast_node_id(node);
  dot << "  \"" << node_id << "\" [shape=" << node->DotShape()
      << ", label=\"" << escape_dot_label(node->DotLabel())
      << "\"];\n";

  const auto children = node->DotChildren();
  bool has_child = false;
  for (const auto &child : children) {
    if (!child) {
      continue;
    }
    has_child = true;
    const std::string child_id = ast_node_id(child);
    dot << "  \"" << node_id << "\" -> \"" << child_id << "\";\n";
    append_ast_subgraph(dot, child, visited, leaf_node_ids);
  }

  if (!has_child) {
    leaf_node_ids.push_back(node_id);
  }
}

void append_leaf_rank_constraints(std::ostringstream &dot,
                                  const std::vector<std::string> &leaf_node_ids) {
  if (leaf_node_ids.empty()) {
    return;
  }

  dot << "  // Keep AST leaves aligned.\n"
      << "  {\n"
      << "    rank=sink;\n";

  if (leaf_node_ids.size() == 1) {
    dot << "    \"" << leaf_node_ids.front() << "\";\n";
  } else {
    dot << "    edge [style=invis, weight=10];\n    ";
    for (std::size_t i = 0; i < leaf_node_ids.size(); ++i) {
      if (i != 0) {
        dot << " -> ";
      }
      dot << "\"" << leaf_node_ids[i] << "\"";
    }
    dot << ";\n";
  }

  dot << "  }\n";
}

} // namespace

export namespace fft_compiler {

class AlgebraicDFTKernelGraphDotWriter {
public:
  [[nodiscard]] std::string ToDot(const AlgebraicDFTKernelGraph &graph) const {
    std::ostringstream dot;
    dot << "digraph algebraic_dft_kernel_graph {\n"
        << "  rankdir=RL;\n"
        << "  graph [labelloc=t, label=\"dft_size=" << graph.dft_size
        << ", dfts_per_kernel=" << graph.dfts_per_kernel << "\"];\n"
        << "  \"graph_root\" [shape=box, label=\"AlgebraicDFTKernelGraph\"];\n";

    for (std::size_t i = 0; i < graph.prime_factor_exponents.size(); ++i) {
      const auto &[prime_factor, exponent] = graph.prime_factor_exponents[i];
      dot << "  \"prime_factor_" << i << "\" [shape=box, label=\"prime "
          << prime_factor << " exponent " << exponent << "\"];\n"
          << "  \"graph_root\" -> \"prime_factor_" << i << "\";\n";
    }

    std::vector<const AlgebraicASTNode *> visited;
    std::vector<std::string> leaf_node_ids;
    for (std::size_t i = 0; i < graph.nodes.size(); ++i) {
      const auto &node = graph.nodes[i];
      const std::string graph_node_id = "graph_node_" + std::to_string(i);
      dot << "  \"" << graph_node_id << "\" [shape=ellipse, label=\"graph.nodes["
          << i << "]\"];\n"
          << "  \"graph_root\" -> \"" << graph_node_id << "\";\n";

      if (!node) {
        dot << "  \"" << graph_node_id << "_null\" [shape=point];\n"
            << "  \"" << graph_node_id << "\" -> \"" << graph_node_id
            << "_null\";\n";
        continue;
      }

      const std::string ast_id = ast_node_id(node);
      dot << "  \"" << graph_node_id << "\" -> \"" << ast_id << "\";\n";
      append_ast_subgraph(dot, node, visited, leaf_node_ids);
    }

    append_leaf_rank_constraints(dot, leaf_node_ids);

    dot << "}\n";
    return dot.str();
  }

  bool WriteDotFile(const AlgebraicDFTKernelGraph &graph,
                    std::string_view output_path) const {
    const std::filesystem::path path(output_path);
    if (path.has_parent_path()) {
      std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
      return false;
    }

    out << ToDot(graph);
    return static_cast<bool>(out);
  }
};

} // namespace fft_compiler
