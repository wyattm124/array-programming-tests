import <charconv>;
import <iostream>;
import <optional>;
import <string_view>;
import <vector>;

import fft.compiler.dot;
import fft.compiler.graph;

namespace {

struct CompilerOptions {
  unsigned int dft_size = 0;
  unsigned int dfts_per_kernel = 1;
  const char *dot_output_path = "build/fft_compiler_graph.dot";
};

std::optional<unsigned int> parse_unsigned(std::string_view arg) {
  unsigned int value = 0;
  const char *begin = arg.data();
  const char *end = begin + arg.size();
  auto [ptr, ec] = std::from_chars(begin, end, value);
  if (ec != std::errc() || ptr != end || value == 0) {
    return std::nullopt;
  }
  return value;
}

void print_usage(const char *program_name) {
  std::cerr << "Usage: " << program_name << " --dft-size <positive integer> "
            << "[--dfts-per-kernel <positive integer>] "
            << "[--dot-out <path>]\n"
            << "  --dft-size is required.\n"
            << "  --dfts-per-kernel is optional and defaults to 1.\n"
            << "  --dot-out is optional and defaults to "
               "build/fft_compiler_graph.dot.\n";
}

std::optional<CompilerOptions> parse_args(int argc, char **argv) {
  CompilerOptions options;
  bool saw_dft_size = false;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];

    auto parse_flag_value = [&](unsigned int &out_value) -> bool {
      if (i + 1 >= argc) {
        return false;
      }
      const auto parsed = parse_unsigned(argv[++i]);
      if (!parsed) {
        return false;
      }
      out_value = *parsed;
      return true;
    };

    if (arg == "--dft-size") {
      if (!parse_flag_value(options.dft_size)) {
        return std::nullopt;
      }
      saw_dft_size = true;
      continue;
    }

    if (arg == "--dfts-per-kernel") {
      if (!parse_flag_value(options.dfts_per_kernel)) {
        return std::nullopt;
      }
      continue;
    }

    if (arg == "--dot-out") {
      if (i + 1 >= argc) {
        return std::nullopt;
      }
      options.dot_output_path = argv[++i];
      continue;
    }

    return std::nullopt;
  }

  if (!saw_dft_size) {
    return std::nullopt;
  }

  return options;
}

} // namespace

int main(int argc, char **argv) {
  const auto options = parse_args(argc, argv);
  if (!options) {
    print_usage(argv[0]);
    return 1;
  }

  const fft_compiler::AlgebraicDFTKernelGraphBuilder graph_builder;
  const fft_compiler::AlgebraicDFTKernelGraph graph =
      graph_builder.Build(options->dft_size, options->dfts_per_kernel);

  const fft_compiler::AlgebraicDFTKernelGraphDotWriter dot_writer;
  if (!dot_writer.WriteDotFile(graph, options->dot_output_path)) {
    std::cerr << "failed to write dot file: " << options->dot_output_path
              << '\n';
    return 1;
  }

  std::cout << "dft_size=" << graph.dft_size << '\n'
            << "dfts_per_kernel=" << graph.dfts_per_kernel << '\n'
            << "dot_file=" << options->dot_output_path << '\n'
            << "prime_factor_exponents=";

  for (std::size_t i = 0; i < graph.prime_factor_exponents.size(); ++i) {
    const auto &factor = graph.prime_factor_exponents[i];
    std::cout << '(' << factor.first << ", " << factor.second << ')';
    if (i + 1 < graph.prime_factor_exponents.size()) {
      std::cout << ' ';
    }
  }
  std::cout << '\n';
  return 0;
}
