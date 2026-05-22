#include "fft_comp_unit.hpp"
#include "../src/arch_config.hpp"
#include "../src/fft.hpp"
#include <benchmark/benchmark.h>
#include <fftw3.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

constexpr const char *kDefaultBenchmarkSpecPath = "bench/specs/fft_bench.yaml";

struct BenchmarkBaseline {
    double cpu_time_ns = 0.0;
};

using BaselineMap = std::unordered_map<std::string, BenchmarkBaseline>;

struct BenchmarkGroup {
    std::string size_label;
    std::string section_label;
};

template <unsigned long N>
void init_case(float (&time_domain)[2 * N], float (&freq_domain)[2 * N]) {
    std::fill_n(time_domain, 2 * N, 0.0f);
    std::fill_n(freq_domain, 2 * N, 0.0f);
    FFT::wave_gen_lcg(time_domain, freq_domain, N);
}

template <unsigned long N>
void interleaved_to_split(const float *interleaved, float *real, float *imag) {
    for (unsigned long i = 0; i < N; i++) {
        real[i] = interleaved[2 * i];
        imag[i] = interleaved[2 * i + 1];
    }
}

double adjusted_cpu_time_ns(const benchmark::BenchmarkReporter::Run &run) {
    const double adjusted = run.GetAdjustedCPUTime();
    switch (run.time_unit) {
        case benchmark::kNanosecond:
            return adjusted;
        case benchmark::kMicrosecond:
            return adjusted * 1e3;
        case benchmark::kMillisecond:
            return adjusted * 1e6;
        case benchmark::kSecond:
            return adjusted * 1e9;
    }
    return adjusted;
}

std::string benchmark_name(unsigned long size, const char *kind,
                           std::optional<std::string> simd = std::nullopt) {
    std::string name = "N" + std::to_string(size) + "/" + kind;
    if (simd.has_value() && *simd != "NONE") {
        name += " [" + *simd + "]";
    }
    return name;
}

std::optional<double> parse_baseline_cpu_time_ns(const YAML::Node &config) {
    if (!config || config.IsNull()) {
        return std::nullopt;
    }
    if (config.IsScalar()) {
        return config.as<double>();
    }
    if (!config.IsMap()) {
        throw std::runtime_error(
            "benchmark config entries must be scalars or mappings.");
    }
    const YAML::Node baseline = config["baseline_cpu_time_ns"];
    if (!baseline) {
        return std::nullopt;
    }
    return baseline.as<double>();
}

template <typename Callback>
void for_each_enabled_variant(const YAML::Node &node, const char *field_name,
                              Callback &&callback) {
    if (!node) {
        return;
    }
    if (node.IsSequence()) {
        for (const auto &entry : node) {
            callback(entry.as<std::string>(), YAML::Node());
        }
        return;
    }
    if (!node.IsMap()) {
        throw std::runtime_error(std::string("'") + field_name +
                                 "' must be a sequence or mapping.");
    }
    for (const auto &entry : node) {
        callback(entry.first.as<std::string>(), entry.second);
    }
}

BaselineMap load_baselines(const std::string &path) {
    BaselineMap baselines;

    std::ifstream stream(path);
    if (!stream.good()) {
        std::cerr << "warning: benchmark spec file not found: " << path
                  << "\n";
        return baselines;
    }

    const YAML::Node root = YAML::Load(stream);
    const YAML::Node benchmarks = root["benchmarks"];
    if (!benchmarks) {
        return baselines;
    }
    if (!benchmarks.IsSequence()) {
        throw std::runtime_error("'benchmarks' must be a YAML sequence.");
    }

    for (const auto &case_node : benchmarks) {
        const unsigned long size = case_node["size"].as<unsigned long>();

        if (const auto comp_unit_baseline =
                parse_baseline_cpu_time_ns(case_node["comp_unit_raw"])) {
            baselines.emplace(benchmark_name(size, "CompUnitRaw"),
                              BenchmarkBaseline{*comp_unit_baseline});
        }

        for_each_enabled_variant(
            case_node["raw_fft"], "raw_fft",
            [&](const std::string &simd, const YAML::Node &config) {
                const auto baseline = parse_baseline_cpu_time_ns(config);
                if (!baseline) {
                    return;
                }
                baselines.emplace(benchmark_name(size, "RawFFT", simd),
                                  BenchmarkBaseline{*baseline});
            });

        if (const auto fftw_baseline =
                parse_baseline_cpu_time_ns(case_node["fftw"])) {
            baselines.emplace(benchmark_name(size, "FFTW"),
                              BenchmarkBaseline{*fftw_baseline});
        }
    }

    return baselines;
}

class BaselineDeltaReporter : public benchmark::ConsoleReporter {
public:
    explicit BaselineDeltaReporter(BaselineMap baselines)
        : benchmark::ConsoleReporter(), baselines_(std::move(baselines)) {}

protected:
    void PrintHeader(const Run &run) override {
        std::ostream &out = GetOutputStream();
        name_field_width_ = std::max<std::size_t>(
            name_field_width_, std::max<std::size_t>(run.benchmark_name().size(), 9));
        const std::size_t total_width = name_field_width_ + kTimeWidth + kCpuWidth +
                                        kIterationsWidth + kBaselineWidth +
                                        kDeltaWidth + kChangeWidth;

        out << std::string(total_width, '-') << '\n'
            << std::left << std::setw(static_cast<int>(name_field_width_))
            << "Benchmark"
            << std::right << std::setw(kTimeWidth) << "Time"
            << std::setw(kCpuWidth) << "CPU"
            << std::setw(kIterationsWidth) << "Iterations"
            << std::setw(kBaselineWidth) << "Baseline"
            << std::setw(kDeltaWidth) << "Delta"
            << std::setw(kChangeWidth) << "Change"
            << '\n'
            << std::string(total_width, '-') << '\n';
    }

    void PrintRunData(const Run &run) override {
        std::ostream &out = GetOutputStream();
        const BenchmarkGroup group = parse_group(run);
        emit_group_headers(group);
        const std::string row_label = display_name(run, group);

        out << std::left << std::setw(static_cast<int>(name_field_width_))
            << row_label
            << colorized_time_cell(format_time(run.GetAdjustedRealTime(), run.time_unit),
                                   kTimeWidth)
            << colorized_time_cell(format_time(run.GetAdjustedCPUTime(), run.time_unit),
                                   kCpuWidth)
            << colorized_info_cell(std::to_string(run.iterations), kIterationsWidth);

        const auto it = baselines_.find(run.benchmark_name());
        if (run.run_type != Run::RT_Iteration || it == baselines_.end() ||
            it->second.cpu_time_ns == 0.0) {
            out << padded_cell("", kBaselineWidth)
                << padded_cell("", kDeltaWidth)
                << padded_cell("", kChangeWidth)
                << '\n';
            return;
        }

        const double current_ns = adjusted_cpu_time_ns(run);
        const double baseline_ns = it->second.cpu_time_ns;
        const double delta_ns = current_ns - baseline_ns;
        const double delta_pct = (delta_ns / baseline_ns) * 100.0;

        out << colorized_info_cell(format_ns(baseline_ns), kBaselineWidth)
            << colorized_cell(format_delta_ns(delta_ns), kDeltaWidth, delta_ns)
            << colorized_cell(format_delta_pct(delta_pct), kChangeWidth, delta_ns)
            << '\n';
    }

private:
    static constexpr int kTimeWidth = 14;
    static constexpr int kCpuWidth = 14;
    static constexpr int kIterationsWidth = 14;
    static constexpr int kBaselineWidth = 14;
    static constexpr int kDeltaWidth = 14;
    static constexpr int kChangeWidth = 12;

    static const char *time_unit_suffix(benchmark::TimeUnit unit) {
        switch (unit) {
            case benchmark::kNanosecond:
                return "ns";
            case benchmark::kMicrosecond:
                return "us";
            case benchmark::kMillisecond:
                return "ms";
            case benchmark::kSecond:
                return "s";
        }
        return "";
    }

    static std::string format_fixed(double value, int precision,
                                    const char *suffix, bool show_pos = false) {
        std::ostringstream stream;
        if (show_pos) {
            stream << std::showpos;
        }
        stream << std::fixed << std::setprecision(precision) << value;
        if (suffix[0] != '\0') {
            stream << ' ' << suffix;
        }
        return stream.str();
    }

    static std::string format_time(double value, benchmark::TimeUnit unit) {
        return format_fixed(value, 2, time_unit_suffix(unit));
    }

    static std::string format_ns(double value) {
        return format_fixed(value, 2, "ns");
    }

    static std::string format_delta_ns(double value) {
        return format_fixed(value, 2, "ns", true);
    }

    static std::string format_delta_pct(double value) {
        return format_fixed(value, 2, "%", true);
    }

    static std::string padded_cell(const std::string &text, int width) {
        std::ostringstream stream;
        stream << std::right;
        stream << std::setw(width) << text;
        return stream.str();
    }

    static BenchmarkGroup parse_group(const Run &run) {
        const std::string function_name =
            run.run_name.function_name.empty() ? run.benchmark_name()
                                               : run.run_name.function_name;
        const std::size_t slash = function_name.find('/');
        if (slash == std::string::npos) {
            return {"Benchmark", function_name};
        }

        std::string size_label = function_name.substr(0, slash);
        if (!size_label.empty() && size_label.front() == 'N') {
            size_label = "FFT Size " + size_label.substr(1);
        }

        return {size_label, function_name.substr(slash + 1)};
    }

    static std::string display_name(const Run &run, const BenchmarkGroup &group) {
        if (run.run_type == Run::RT_Iteration) {
            return "  " + group.section_label;
        }
        if (run.run_type == Run::RT_Aggregate && !run.aggregate_name.empty()) {
            return "    " + run.aggregate_name;
        }
        return "";
    }

    void emit_group_headers(const BenchmarkGroup &group) {
        std::ostream &out = GetOutputStream();

        if (group.size_label != current_size_label_) {
            if (!current_size_label_.empty()) {
                out << '\n';
            }
            out << group.size_label << '\n';
            current_size_label_ = group.size_label;
        }
    }

    std::string colorized_time_cell(const std::string &text, int width) const {
        const std::string padded = padded_cell(text, width);
        if (!(output_options_ & benchmark::ConsoleReporter::OO_Color)) {
            return padded;
        }
        return std::string("\033[33m") + padded + "\033[0m";
    }

    std::string colorized_info_cell(const std::string &text, int width) const {
        const std::string padded = padded_cell(text, width);
        if (!(output_options_ & benchmark::ConsoleReporter::OO_Color)) {
            return padded;
        }
        return std::string("\033[36m") + padded + "\033[0m";
    }

    std::string colorized_cell(const std::string &text, int width,
                               double delta_ns) const {
        const std::string padded = padded_cell(text, width);
        if (!(output_options_ & benchmark::ConsoleReporter::OO_Color) ||
            delta_ns == 0.0) {
            return padded;
        }

        const char *color = delta_ns < 0.0 ? "\033[32m" : "\033[31m";
        return std::string(color) + padded + "\033[0m";
    }

    BaselineMap baselines_;
    std::string current_size_label_;
};

template <unsigned long N>
void process_fftw(benchmark::State &state, const float *time_domain) {
    fftwf_iodim dims[1] = {{static_cast<int>(N), 1, 1}};
    float *in_real = static_cast<float *>(fftwf_malloc(sizeof(float) * N));
    float *in_imag = static_cast<float *>(fftwf_malloc(sizeof(float) * N));
    float *out_real = static_cast<float *>(fftwf_malloc(sizeof(float) * N));
    float *out_imag = static_cast<float *>(fftwf_malloc(sizeof(float) * N));
    interleaved_to_split<N>(time_domain, in_real, in_imag);
    fftwf_plan p = fftwf_plan_guru_split_dft(
        1, dims, 0, nullptr, in_real, in_imag, out_real, out_imag, FFTW_MEASURE);

    for (auto _ : state) {
        fftwf_execute(p);
    }

    fftwf_destroy_plan(p);
    fftwf_free(in_real);
    fftwf_free(in_imag);
    fftwf_free(out_real);
    fftwf_free(out_imag);
}

template <unsigned long N, SIMD_TYPE simd>
void process_fft(benchmark::State &state, const float *time_domain) {
    alignas(MY_MAX_ALIGNMENT) float temp_time_domain[2 * N];
    alignas(MY_MAX_ALIGNMENT) float temp_freq_domain[2 * N] = {0};
    std::copy_n(time_domain, 2 * N, temp_time_domain);
    FFT::FFTPlan<float, simd>::template Init<N>();

    for (auto _ : state) {
        FFT::FFTPlan<float, simd>::template fft<N>(temp_time_domain, temp_freq_domain);
    }
}

template <unsigned long N>
void process_comp_unit(benchmark::State &state, const float *time_domain) {
    static_assert(N >= 2 && N <= 8, "Comp unit benchmark only supports sizes 2 through 8.");
    alignas(MY_MAX_ALIGNMENT) float temp_time_domain[2 * N];
    alignas(MY_MAX_ALIGNMENT) float temp_freq_domain[2 * N] = {0};
    std::copy_n(time_domain, 2 * N, temp_time_domain);

    for (auto _ : state) {
        if constexpr (N == 2) {
            FFT::fft_2(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 3) {
            FFT::fft_3(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 4) {
            FFT::fft_4(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 5) {
            FFT::fft_5(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 6) {
            FFT::fft_6(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 7) {
            FFT::fft_7(temp_time_domain, temp_freq_domain);
        } else if constexpr (N == 8) {
            FFT::fft_8(temp_time_domain, temp_freq_domain);
        }
    }
}

template <unsigned long N>
void register_comp_unit_benchmark(const char *name) {
    benchmark::RegisterBenchmark(name, [](benchmark::State &state) {
        alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N];
        alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N];
        init_case<N>(time_domain, freq_domain);
        process_comp_unit<N>(state, time_domain);
    });
}

template <unsigned long N, SIMD_TYPE simd>
void register_fft_benchmark(const char *name) {
    benchmark::RegisterBenchmark(name, [](benchmark::State &state) {
        alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N];
        alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N];
        init_case<N>(time_domain, freq_domain);
        process_fft<N, simd>(state, time_domain);
    });
}

template <unsigned long N>
void register_fftw_benchmark(const char *name) {
    benchmark::RegisterBenchmark(name, [](benchmark::State &state) {
        alignas(MY_MAX_ALIGNMENT) float time_domain[2 * N];
        alignas(MY_MAX_ALIGNMENT) float freq_domain[2 * N];
        init_case<N>(time_domain, freq_domain);
        process_fftw<N>(state, time_domain);
    });
}

void register_generated_benchmarks() {
#include "generated/fft_bench_none.inc"
#if FFT_HAS_AVX2
#include "generated/fft_bench_avx2.inc"
#endif
#if FFT_HAS_NEON
#include "generated/fft_bench_neon.inc"
#endif
}

} // namespace

int main(int argc, char **argv) {
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }

    register_generated_benchmarks();

    const char *spec_path_env = std::getenv("FFT_BENCH_SPEC");
    const std::string spec_path =
        spec_path_env != nullptr ? spec_path_env : kDefaultBenchmarkSpecPath;

    BaselineMap baselines;
    try {
        baselines = load_baselines(spec_path);
    } catch (const std::exception &ex) {
        std::cerr << "error: failed to load benchmark data from "
                  << spec_path << ": " << ex.what() << "\n";
        return 1;
    }

    BaselineDeltaReporter reporter(std::move(baselines));
    benchmark::RunSpecifiedBenchmarks(&reporter);
    benchmark::Shutdown();
    return 0;
}
