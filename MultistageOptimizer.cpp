/*
Filename: MultistageOptimizer.cpp
Description : Calculate optimal mass for each stage given a DeltaV requirement on a multistage Rocket
Source of math equations : http://www.projectrho.com/public_html/rocket/multistage.php
Author: Florian Zimmer
Version : 0.8
Last Changed : 31.03.2021
Changelog :
0.1 basic functionality/port from python project
0.2 debugging
0.3 json input
0.4 cleanup code
0.5 fix Perfomance Issues and Multithreading
0.6 split files
0.7 added config
0.8 added icon
*/

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#endif

#include "rocket_definition/Stage.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "math/math.hpp"
#include "math/massCalc.hpp"
#include "math/validation.hpp"
#include "optimizer/budget.hpp"
#include "optimizer/solver.hpp"
#include "optimizer/reporting.hpp"
#include "read_json/read_json.hpp"
#include "util/config_parsing.hpp"
#include "util/format_numbers.hpp"
#include "util/text_table.hpp"

using namespace std;
using MassType = double;

struct MemoryCounters {
    std::uint64_t working_set_bytes = 0;
    std::uint64_t private_bytes = 0;
};

static bool readCurrentProcessMemory(MemoryCounters& out)
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    std::memset(&pmc, 0, sizeof(pmc));
    if (!K32GetProcessMemoryInfo(GetCurrentProcess(),
                                reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
                                static_cast<DWORD>(sizeof(pmc)))) {
        return false;
    }
    out.working_set_bytes = static_cast<std::uint64_t>(pmc.WorkingSetSize);
    out.private_bytes = static_cast<std::uint64_t>(pmc.PrivateUsage);
    return true;
#else
    (void)out;
    return false;
#endif
}

template <typename Fn>
static bool runWithPeakMemory(Fn&& fn, MemoryCounters& peak_out, int sample_ms = 5)
{
    std::atomic<bool> stop{false};
    MemoryCounters peak;

    if (sample_ms < 0) {
        sample_ms = 0;
    }

    std::thread sampler([&] {
        for (;;) {
            MemoryCounters current;
            if (readCurrentProcessMemory(current)) {
                peak.working_set_bytes = std::max(peak.working_set_bytes, current.working_set_bytes);
                peak.private_bytes = std::max(peak.private_bytes, current.private_bytes);
            }
            if (stop.load(std::memory_order_relaxed)) {
                break;
            }
            if (sample_ms > 0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(sample_ms));
            }
        }
    });

    bool ok = false;
    try {
        ok = fn();
    }
    catch (...) {
        stop.store(true, std::memory_order_relaxed);
        sampler.join();
        throw;
    }

    stop.store(true, std::memory_order_relaxed);
    sampler.join();
    peak_out = peak;
    return ok;
}

struct BenchmarkTimings {
    double fill_seconds = 0.0;
    double tuple_seconds = 0.0;
    double mass_seconds = 0.0;
    double min_seconds = 0.0;
    double total_seconds = 0.0;
};

struct OptimizationSummary {
    BenchmarkTimings timings;
    unsigned long long n_combinations = 0;
    MassType min_mass = 0.0;
    long long min_index = 0;
    int stage_count = 0;
    int precision = 0;
    int threads_used = 0;
    int deltaVmission = 0;
    std::vector<unsigned char> best_distribution;
};

struct RunOptions {
    bool print_output = true;
    bool pause_on_exit = true;
    bool force_verbose_off = false;
    int threads_override = 0;
    string full_grid_mode_override;
};

struct ProgramOptions {
    bool benchmark = false;
    bool benchmark_compare = false;
    int benchmark_iterations = 1;
    int benchmark_warmup = 0;
    int benchmark_threads = 0;
    double benchmark_max_seconds = -1.0;
    string benchmark_compare_format = "json";
    string benchmark_csv_path = "benchmark_results.csv";
    string full_grid_mode_override;
    bool no_prompt = false;
    bool show_help = false;
    string config_path = "config/config.json";
    string benchmark_config_path;
};

using Clock = std::chrono::steady_clock;

static double elapsedSeconds(Clock::time_point start, Clock::time_point end)
{
    return std::chrono::duration<double>(end - start).count();
}

static std::string formatBytesHuman(std::uint64_t bytes)
{
    std::ostringstream oss;
    constexpr double gib_div = 1024.0 * 1024.0 * 1024.0;
    double gib = static_cast<double>(bytes) / gib_div;
    oss << std::fixed << std::setprecision(2) << gib << " GiB (" << bytes << " bytes)";
    return oss.str();
}

static void output2Darr(unsigned char* arr, int sizeY, int sizeX)
{
    for (int i = 0; i < sizeY; i++) {
        for (int j = 0; j < sizeX; j++) {
            std::cout << (int)arr[j * sizeY + i] << "\t";
        }
        std::cout << "\n";
    }
}

static void outputArr(double* arr, int sizeX)
{
    for (int i = 0; i < sizeX; i++) {
        std::cout << arr[i] << "%\t";
    }
}

static string trimLine(const string& value)
{
    size_t start = 0;
    while (start < value.size() && (value[start] == ' ' || value[start] == '\t' || value[start] == '\r' || value[start] == '\n')) {
        start++;
    }
    size_t end = value.size();
    while (end > start && (value[end - 1] == ' ' || value[end - 1] == '\t' || value[end - 1] == '\r' || value[end - 1] == '\n')) {
        end--;
    }
    return value.substr(start, end - start);
}

static bool readFirstLine(const string& path, string& line)
{
    ifstream file(path);
    if (!file.is_open()) {
        return false;
    }
    if (!std::getline(file, line)) {
        return false;
    }
    line = trimLine(line);
    return !line.empty();
}

static string readPackedRef(const string& packed_refs_path, const string& ref)
{
    ifstream file(packed_refs_path);
    if (!file.is_open()) {
        return "";
    }
    string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '^') {
            continue;
        }
        std::istringstream iss(line);
        string hash;
        string name;
        if (!(iss >> hash >> name)) {
            continue;
        }
        if (name == ref) {
            return trimLine(hash);
        }
    }
    return "";
}

static bool readEnvVar(const char* name, string& value)
{
#if defined(_MSC_VER)
    char* buffer = nullptr;
    size_t length = 0;
    if (_dupenv_s(&buffer, &length, name) != 0 || buffer == nullptr) {
        return false;
    }
    value.assign(buffer);
    free(buffer);
    value = trimLine(value);
    return !value.empty();
#else
    const char* result = std::getenv(name);
    if (!result || !*result) {
        return false;
    }
    value = trimLine(result);
    return !value.empty();
#endif
}

static string getGitHead()
{
    string value;
    if (readEnvVar("GIT_HEAD", value)) {
        return value;
    }
    if (readEnvVar("GITHUB_SHA", value)) {
        return value;
    }
    if (readEnvVar("CI_COMMIT_SHA", value)) {
        return value;
    }

    string head_line;
    if (!readFirstLine(".git/HEAD", head_line)) {
        return "unknown";
    }
    const string ref_prefix = "ref:";
    if (head_line.compare(0, ref_prefix.size(), ref_prefix) == 0) {
        string ref_path = trimLine(head_line.substr(ref_prefix.size()));
        if (!ref_path.empty() && ref_path[0] == ':') {
            ref_path = trimLine(ref_path.substr(1));
        }
        string ref_line;
        if (readFirstLine(".git/" + ref_path, ref_line)) {
            return ref_line;
        }
        string packed_ref = readPackedRef(".git/packed-refs", ref_path);
        if (!packed_ref.empty()) {
            return packed_ref;
        }
        return "unknown";
    }
    return head_line;
}

static string jsonEscape(const string& value)
{
    string out;
    out.reserve(value.size());
    for (char c : value) {
        switch (c) {
        case '\\':
            out += "\\\\";
            break;
        case '"':
            out += "\\\"";
            break;
        case '\b':
            out += "\\b";
            break;
        case '\f':
            out += "\\f";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (static_cast<unsigned char>(c) < 0x20) {
                std::ostringstream oss;
                oss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)static_cast<unsigned char>(c);
                out += oss.str();
            }
            else {
                out += c;
            }
            break;
        }
    }
    return out;
}

static bool parseInt(const char* value, int& out)
{
    if (!value || *value == '\0') {
        return false;
    }
    char* end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (*end != '\0') {
        return false;
    }
    if (parsed > INT_MAX || parsed < INT_MIN) {
        return false;
    }
    out = static_cast<int>(parsed);
    return true;
}

static bool parseDouble(const char* value, double& out)
{
    if (!value || *value == '\0') {
        return false;
    }
    char* end = nullptr;
    double parsed = std::strtod(value, &end);
    if (*end != '\0') {
        return false;
    }
    out = parsed;
    return true;
}

static bool parseUInt64Json(const json& value, std::uint64_t& out)
{
    if (value.is_number_unsigned()) {
        out = value.get<std::uint64_t>();
        return true;
    }
    if (value.is_number_integer()) {
        long long signed_value = value.get<long long>();
        if (signed_value < 0) {
            return false;
        }
        out = static_cast<std::uint64_t>(signed_value);
        return true;
    }
    if (value.is_number_float()) {
        double number = value.get<double>();
        if (!std::isfinite(number) || number < 0.0 || number > static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
            return false;
        }
        out = static_cast<std::uint64_t>(number);
        return true;
    }
    return false;
}

static void printUsage()
{
    std::cout << "Usage: MultistageOptimizer.exe [options]\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --config <path>                 Config file for normal runs.\n";
    std::cout << "  --no-prompt                      Skip waiting for Enter before exit.\n";
    std::cout << "  --benchmark                      Run performance benchmark.\n";
    std::cout << "  --benchmark-compare              Run multiple benchmark strategies and compare.\n";
    std::cout << "  --benchmark-compare-format <f>   Output format for --benchmark-compare (json|table).\n";
    std::cout << "  --benchmark-config <path>        Config file for benchmark runs.\n";
    std::cout << "  --benchmark-iterations <count>   Number of benchmark samples.\n";
    std::cout << "  --benchmark-warmup <count>       Warmup runs before sampling.\n";
    std::cout << "  --benchmark-threads <count>      Override thread count.\n";
    std::cout << "  --fullgrid-mode <mode>           Full-grid mode: materialize|streaming (default: materialize).\n";
    std::cout << "  --benchmark-max-seconds <sec>    Fail if average total exceeds this.\n";
    std::cout << "  --benchmark-csv <path>           Append benchmark results to CSV.\n";
    std::cout << "  --help                           Show this help.\n";
}

static bool parseOptions(int argc, char** argv, ProgramOptions& options)
{
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            options.show_help = true;
            return true;
        }
        if (arg == "--config") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --config.\n";
                return false;
            }
            options.config_path = argv[++i];
        }
        else if (arg == "--no-prompt") {
            options.no_prompt = true;
        }
        else if (arg == "--benchmark") {
            options.benchmark = true;
        }
        else if (arg == "--benchmark-compare") {
            options.benchmark_compare = true;
        }
        else if (arg == "--benchmark-compare-format") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-compare-format.\n";
                return false;
            }
            options.benchmark_compare_format = argv[++i];
            if (options.benchmark_compare_format != "json" && options.benchmark_compare_format != "table") {
                std::cerr << "Invalid --benchmark-compare-format value (expected json|table).\n";
                return false;
            }
        }
        else if (arg == "--benchmark-config") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-config.\n";
                return false;
            }
            options.benchmark_config_path = argv[++i];
        }
        else if (arg == "--benchmark-iterations") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-iterations.\n";
                return false;
            }
            int value = 0;
            if (!parseInt(argv[++i], value) || value < 1) {
                std::cerr << "Invalid --benchmark-iterations value.\n";
                return false;
            }
            options.benchmark_iterations = value;
        }
        else if (arg == "--benchmark-warmup") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-warmup.\n";
                return false;
            }
            int value = 0;
            if (!parseInt(argv[++i], value) || value < 0) {
                std::cerr << "Invalid --benchmark-warmup value.\n";
                return false;
            }
            options.benchmark_warmup = value;
        }
        else if (arg == "--benchmark-threads") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-threads.\n";
                return false;
            }
            int value = 0;
            if (!parseInt(argv[++i], value) || value < 1) {
                std::cerr << "Invalid --benchmark-threads value.\n";
                return false;
            }
            options.benchmark_threads = value;
        }
        else if (arg == "--fullgrid-mode") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --fullgrid-mode.\n";
                return false;
            }
            options.full_grid_mode_override = argv[++i];
        }
        else if (arg == "--benchmark-max-seconds") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-max-seconds.\n";
                return false;
            }
            double value = 0.0;
            if (!parseDouble(argv[++i], value) || value <= 0.0) {
                std::cerr << "Invalid --benchmark-max-seconds value.\n";
                return false;
            }
            options.benchmark_max_seconds = value;
        }
        else if (arg == "--benchmark-csv") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for --benchmark-csv.\n";
                return false;
            }
            options.benchmark_csv_path = argv[++i];
        }
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return false;
        }
    }
    return true;
}

static int resolveThreadCount(bool useMultiCore, int threads_override)
{
    if (threads_override > 0) {
        return threads_override;
    }
    if (!useMultiCore) {
        return 1;
    }
    unsigned int hw_threads = std::thread::hardware_concurrency();
    if (hw_threads == 0) {
        return 1;
    }
    return static_cast<int>(hw_threads);
}

static bool parseFullGridModeString(const std::string& value, optimizer::FullGridMode& out)
{
    std::string s;
    s.reserve(value.size());
    for (char c : value) {
        s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }

    if (s == "materialize" || s == "materialized") {
        out = optimizer::FullGridMode::Materialize;
        return true;
    }
    if (s == "streaming" || s == "stream") {
        out = optimizer::FullGridMode::Streaming;
        return true;
    }

    return false;
}

static bool runOptimization(const json& config, const RunOptions& options, OptimizationSummary& summary, std::string* error_out = nullptr)
{
    if (error_out) {
        error_out->clear();
    }
    auto setError = [&](const std::string& message) {
        if (error_out) {
            *error_out = message;
        }
    };

    bool verbose = false;
    if (config.contains("verbose")) {
        if (!config["verbose"].is_boolean()) {
            std::cerr << "Invalid verbose value in config.\n";
            setError("Invalid verbose value in config.");
            return false;
        }
        verbose = config["verbose"].get<bool>();
    }
    if (options.force_verbose_off) {
        verbose = false;
    }

    optimizer::FullGridMode full_grid_mode = optimizer::FullGridMode::Materialize;
    if (config.contains("fullGridMode")) {
        if (!config["fullGridMode"].is_string()) {
            std::cerr << "Invalid fullGridMode value in config.\n";
            setError("Invalid fullGridMode value in config.");
            return false;
        }
        std::string mode = config["fullGridMode"].get<std::string>();
        if (!parseFullGridModeString(mode, full_grid_mode)) {
            std::cerr << "Invalid fullGridMode value in config (expected materialize|streaming).\n";
            setError("Invalid fullGridMode value in config.");
            return false;
        }
    }
    if (!options.full_grid_mode_override.empty()) {
        if (!parseFullGridModeString(options.full_grid_mode_override, full_grid_mode)) {
            std::cerr << "Invalid --fullgrid-mode value (expected materialize|streaming).\n";
            setError("Invalid --fullgrid-mode value.");
            return false;
        }
    }

    bool useMultiCore = false;
    if (config.contains("useMultiCore")) {
        if (!config["useMultiCore"].is_boolean()) {
            std::cerr << "Invalid useMultiCore value in config.\n";
            setError("Invalid useMultiCore value in config.");
            return false;
        }
        useMultiCore = config["useMultiCore"].get<bool>();
    }
    unsigned long long maxRAM = 0;
    if (!config.contains("maxRAM")) {
        std::cerr << "Missing maxRAM value in config.\n";
        setError("Missing maxRAM value in config.");
        return false;
    }
    if (!parseMaxRam(config["maxRAM"], maxRAM)) {
        std::cerr << "Invalid maxRAM value in config.\n";
        setError("Invalid maxRAM value in config.");
        return false;
    }
    if (!config.contains("precision") || !config["precision"].is_number()) {
        std::cerr << "Missing precision value in config.\n";
        setError("Missing precision value in config.");
        return false;
    }
    int precision = config["precision"].get<int>();
    if (precision <= 0 || precision > 255) {
        std::cerr << "Invalid precision value in config.\n";
        setError("Invalid precision value in config.");
        return false;
    }

    vector<Engine> engineList;
    Rocket rocket;
    try {
        engineList = importEngines(config.at("enginesPath").get<string>());
        rocket = importRocket(config.at("rocketPath").get<string>(), engineList);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        setError(e.what());
        return false;
    }

    const std::size_t stage_count = rocket.stages.size();
    if (stage_count < 1) {
        std::cerr << "Rocket has no stages.\n";
        setError("Rocket has no stages.");
        return false;
    }

    int nThreads = resolveThreadCount(useMultiCore, options.threads_override);
    if (options.threads_override > 0) {
        useMultiCore = (nThreads > 1);
    }

    if (nThreads < 1) {
        nThreads = 1;
    }

    bool zoom_enabled = false;
    bool zoom_has_coarse_precision = false;
    bool zoom_has_fine_precision = false;
    optimizer::ZoomOptions zoom_options;
    if (config.contains("zoom") && config["zoom"].is_object()) {
        const auto& zoom = config["zoom"];
        if (zoom.contains("enabled")) {
            zoom_enabled = zoom["enabled"].get<bool>();
        }
        if (zoom.contains("coarse_precision")) {
            zoom_has_coarse_precision = true;
            zoom_options.coarse_precision = zoom["coarse_precision"].get<int>();
        }
        if (zoom.contains("fine_precision")) {
            zoom_has_fine_precision = true;
            zoom_options.fine_precision = zoom["fine_precision"].get<int>();
        }
        if (zoom.contains("window_coarse_steps")) {
            zoom_options.window_coarse_steps = zoom["window_coarse_steps"].get<int>();
        }
        if (zoom.contains("top_k")) {
            zoom_options.top_k = static_cast<std::size_t>(std::max(1, zoom["top_k"].get<int>()));
        }
    }

    zoom_options.verbose = verbose;
    zoom_options.threads = nThreads;

    double reporting_min_stage_dv_mps = 0.0;
    if (config.contains("reporting") && config["reporting"].is_object()) {
        const auto& reporting = config["reporting"];
        if (reporting.contains("min_stage_dv_mps")) {
            if (!reporting["min_stage_dv_mps"].is_number()) {
                std::cerr << "Invalid reporting.min_stage_dv_mps value in config.\n";
                setError("Invalid reporting.min_stage_dv_mps value in config.");
                return false;
            }
            reporting_min_stage_dv_mps = reporting["min_stage_dv_mps"].get<double>();
            if (!std::isfinite(reporting_min_stage_dv_mps) || reporting_min_stage_dv_mps < 0.0) {
                std::cerr << "Invalid reporting.min_stage_dv_mps value in config.\n";
                setError("Invalid reporting.min_stage_dv_mps value in config.");
                return false;
            }
        }
    }

    std::uint64_t max_combinations_budget = 0;
    if (config.contains("maxCombinations")) {
        std::uint64_t parsed = 0;
        if (!parseUInt64Json(config["maxCombinations"], parsed) || parsed == 0) {
            std::cerr << "Invalid maxCombinations value in config.\n";
            setError("Invalid maxCombinations value in config.");
            return false;
        }
        max_combinations_budget = parsed;
    }

    double target_seconds = -1.0;
    if (max_combinations_budget == 0 && config.contains("targetSeconds")) {
        if (!config["targetSeconds"].is_number()) {
            std::cerr << "Invalid targetSeconds value in config.\n";
            setError("Invalid targetSeconds value in config.");
            return false;
        }
        target_seconds = config["targetSeconds"].get<double>();
        if (!std::isfinite(target_seconds) || target_seconds <= 0.0) {
            std::cerr << "Invalid targetSeconds value in config.\n";
            setError("Invalid targetSeconds value in config.");
            return false;
        }

        int calib_precision = precision;
        if (zoom_enabled && zoom_has_fine_precision && zoom_options.fine_precision > 0) {
            calib_precision = zoom_options.fine_precision;
        }
        calib_precision = std::clamp(calib_precision, static_cast<int>(stage_count), 255);
        std::uint64_t calib_total = nCr(calib_precision - 1, static_cast<int>(stage_count) - 1);
        if (calib_total == 0) {
            std::cerr << "Unable to calibrate runtime (invalid precision/stages).\n";
            setError("Unable to calibrate runtime (invalid precision/stages).");
            return false;
        }
        std::uint64_t sample = std::min<std::uint64_t>(5000, calib_total);

        std::vector<unsigned char> dist(sample * stage_count, 1);
        std::uint64_t gen = createTuple(dist.data(), stage_count, sample, calib_precision);
        if (gen != sample) {
            std::cerr << "Unable to calibrate runtime (tuple generation failed).\n";
            setError("Unable to calibrate runtime (tuple generation failed).");
            return false;
        }
        std::vector<double> masses(sample, 0.0);
        auto calib_start = Clock::now();
        if (nThreads == 1) {
            distMass(false, 0, sample, sample, rocket, masses.data(), dist.data(), calib_precision);
        }
        else {
            std::uint64_t threadBlock = sample / static_cast<std::uint64_t>(nThreads);
            std::vector<std::thread> threads;
            threads.reserve(static_cast<std::size_t>(nThreads));
            for (int t = 0; t < nThreads; t++) {
                std::uint64_t begin = static_cast<std::uint64_t>(t) * threadBlock;
                std::uint64_t end = (t == nThreads - 1) ? sample : begin + threadBlock;
                threads.push_back(std::thread([begin, end, sample, &rocket, &masses, &dist, calib_precision] {
                    distMass(false, begin, end, sample, rocket, masses.data(), dist.data(), calib_precision);
                }));
            }
            for (auto& th : threads) {
                th.join();
            }
        }
        auto calib_end = Clock::now();
        double secs = elapsedSeconds(calib_start, calib_end);
        if (!std::isfinite(secs) || secs <= 0.0) {
            std::cerr << "Unable to calibrate runtime (timing failed).\n";
            setError("Unable to calibrate runtime (timing failed).");
            return false;
        }
        double combos_per_sec = static_cast<double>(sample) / secs;
        std::uint64_t estimated = static_cast<std::uint64_t>(std::floor(target_seconds * combos_per_sec));
        if (estimated == 0) {
            estimated = 1;
        }
        max_combinations_budget = estimated;
        if (options.print_output) {
            std::cout << "Calibration: " << combos_per_sec << " combos/sec (" << sample << " samples, precision " << calib_precision << ", threads " << nThreads << ")\n";
            std::cout << "TargetSeconds: " << target_seconds << " => maxCombinations=" << max_combinations_budget << "\n\n";
        }
    }

    int requested_fine_precision = zoom_options.fine_precision;
    if (!zoom_enabled || requested_fine_precision <= 0) {
        requested_fine_precision = precision;
    }

    int effective_precision = requested_fine_precision;
    if (effective_precision <= 0 || effective_precision > 255) {
        std::cerr << "Invalid precision value.\n";
        setError("Invalid precision value.");
        return false;
    }

    int coarse_precision = effective_precision;
    std::uint64_t refinement_budget = std::numeric_limits<std::uint64_t>::max();

    if (max_combinations_budget > 0) {
        if (!zoom_enabled) {
            int chosen =
                (full_grid_mode == optimizer::FullGridMode::Streaming)
                    ? optimizer::choosePrecisionForMaxCombinationsStreaming(static_cast<int>(stage_count),
                                                                            max_combinations_budget,
                                                                            maxRAM,
                                                                            sizeof(MassType),
                                                                            255)
                    : optimizer::choosePrecisionForMaxCombinations(static_cast<int>(stage_count),
                                                                  max_combinations_budget,
                                                                  maxRAM,
                                                                  sizeof(MassType),
                                                                  255);
            if (chosen <= 0) {
                std::cerr << "Unable to satisfy maxCombinations budget with current stage count and RAM.\n";
                setError("Unable to satisfy maxCombinations budget with current stage count and RAM.");
                return false;
            }
            effective_precision = chosen;
            coarse_precision = chosen;
        }
        else {
            const std::uint64_t coarse_budget = std::max<std::uint64_t>(1, static_cast<std::uint64_t>(std::floor(static_cast<double>(max_combinations_budget) * 0.7)));
            int coarse_cap =
                (full_grid_mode == optimizer::FullGridMode::Streaming)
                    ? optimizer::choosePrecisionForMaxCombinationsStreaming(static_cast<int>(stage_count),
                                                                            coarse_budget,
                                                                            maxRAM,
                                                                            sizeof(MassType),
                                                                            effective_precision)
                    : optimizer::choosePrecisionForMaxCombinations(static_cast<int>(stage_count),
                                                                  coarse_budget,
                                                                  maxRAM,
                                                                  sizeof(MassType),
                                                                  effective_precision);
            if (coarse_cap <= 0) {
                std::cerr << "Unable to satisfy maxCombinations coarse budget with current stage count and RAM.\n";
                setError("Unable to satisfy maxCombinations coarse budget with current stage count and RAM.");
                return false;
            }
            int requested_coarse = zoom_options.coarse_precision;
            if (zoom_has_coarse_precision && requested_coarse > 0) {
                coarse_precision = std::min(requested_coarse, coarse_cap);
            }
            else {
                coarse_precision = coarse_cap;
            }

            coarse_precision = std::max(coarse_precision, static_cast<int>(stage_count));
            coarse_precision = std::min(coarse_precision, effective_precision);

            std::uint64_t coarse_evaluated = nCr(coarse_precision - 1, static_cast<int>(stage_count) - 1);
            refinement_budget = (max_combinations_budget > coarse_evaluated) ? (max_combinations_budget - coarse_evaluated) : 1;
        }
    }
    else if (zoom_enabled) {
        coarse_precision = zoom_options.coarse_precision;
        if (coarse_precision <= 0) {
            coarse_precision = std::max(static_cast<int>(stage_count), effective_precision / 2);
        }
        coarse_precision = std::max(coarse_precision, static_cast<int>(stage_count));
        coarse_precision = std::min(coarse_precision, effective_precision);
    }

    std::uint64_t coarse_combinations = nCr(coarse_precision - 1, static_cast<int>(stage_count) - 1);
    if (coarse_combinations == 0) {
        std::cerr << "No combinations generated.\n";
        setError("No combinations generated.");
        return false;
    }

    unsigned long long coarse_required =
        (full_grid_mode == optimizer::FullGridMode::Streaming)
            ? requiredBytesForPrecisionStreaming(coarse_precision, static_cast<int>(stage_count), sizeof(MassType))
            : requiredBytesForPrecision(coarse_precision, static_cast<int>(stage_count), sizeof(MassType));
    if (!zoom_enabled) {
        if (coarse_required >= maxRAM) {
            int max_prec =
                (full_grid_mode == optimizer::FullGridMode::Streaming)
                    ? calcMaxPrecUpStreaming(coarse_precision, static_cast<int>(stage_count), maxRAM, sizeof(MassType))
                    : calcMaxPrecUp(coarse_precision, static_cast<int>(stage_count), maxRAM, sizeof(MassType));
            std::cerr << "Precision too high for the amount of stages and maxRAM.\n";
            std::cerr << "Requested precision: " << coarse_precision << "\n";
            std::cerr << "Stages: " << stage_count << "\n";
            std::cerr << "maxRAM: " << formatBytesHuman(static_cast<std::uint64_t>(maxRAM)) << "\n";
            std::cerr << "Estimated required RAM (full-grid): " << formatBytesHuman(static_cast<std::uint64_t>(coarse_required)) << "\n";
            std::cerr << "Maximum precision for current stages and RAM: " << max_prec << "\n";

            setError("precision too high (requested=" + std::to_string(coarse_precision) + ", max=" + std::to_string(max_prec) + ")");
            return false;
        }
    }
    else {
        if (coarse_required >= maxRAM) {
            int max_prec =
                (full_grid_mode == optimizer::FullGridMode::Streaming)
                    ? calcMaxPrecUpStreaming(coarse_precision, static_cast<int>(stage_count), maxRAM, sizeof(MassType))
                    : calcMaxPrecUp(coarse_precision, static_cast<int>(stage_count), maxRAM, sizeof(MassType));
            std::cerr << "Zoom coarse_precision too high for the amount of stages and maxRAM.\n";
            std::cerr << "Requested coarse_precision: " << coarse_precision << "\n";
            std::cerr << "Stages: " << stage_count << "\n";
            std::cerr << "maxRAM: " << formatBytesHuman(static_cast<std::uint64_t>(maxRAM)) << "\n";
            std::cerr << "Estimated required RAM (coarse full-grid): " << formatBytesHuman(static_cast<std::uint64_t>(coarse_required)) << "\n";
            std::cerr << "Maximum coarse_precision for current stages and RAM: " << max_prec << "\n";

            setError("coarse_precision too high (requested=" + std::to_string(coarse_precision) + ", max=" + std::to_string(max_prec) + ")");
            return false;
        }
    }

    if (options.print_output) {
        std::cout << "For the number of stages and RAM currently given, the maximum precision can be raised to: "
                  << ((full_grid_mode == optimizer::FullGridMode::Streaming)
                          ? calcMaxPrecDownStreaming(255, static_cast<int>(stage_count), maxRAM, sizeof(MassType))
                          : calcMaxPrecDown(255, static_cast<int>(stage_count), maxRAM, sizeof(MassType)))
                  << "\n\n";
        if (zoom_enabled) {
            std::cout << "Zoom enabled: coarse_precision=" << coarse_precision
                      << " fine_precision=" << effective_precision << "\n";
        }
        if (max_combinations_budget > 0) {
            std::cout << "Budget enabled: maxCombinations=" << max_combinations_budget << "\n";
        }
        std::cout << "# Combinations (coarse): " << coarse_combinations << "\n\n";
    }

    auto total_start = Clock::now();

    MassType min = std::numeric_limits<MassType>::infinity();
    long long index = 0;
    std::vector<unsigned char> bestDistro(stage_count, 0);
    BenchmarkTimings timings;

    if (options.print_output) {
        if (useMultiCore) {
            cout << "using " << nThreads << " threads\n";
        }
        else {
            cout << "using 1 thread\n";
        }
    }

    if (!zoom_enabled) {
        optimizer::FullGridResult solve = optimizer::solveFullGrid(rocket, effective_precision, verbose, nThreads, 1, full_grid_mode);
        min = solve.best_mass;
        index = static_cast<long long>(solve.best_index);
        bestDistro = solve.best_units;
        timings.fill_seconds = solve.timings.fill_seconds;
        timings.tuple_seconds = solve.timings.tuple_seconds;
        timings.mass_seconds = solve.timings.mass_seconds;
        timings.min_seconds = solve.timings.min_seconds;
        summary.n_combinations = solve.combinations;
    }
    else {
        zoom_options.coarse_precision = coarse_precision;
        zoom_options.fine_precision = effective_precision;
        zoom_options.max_evaluations = refinement_budget;
        zoom_options.full_grid_mode = full_grid_mode;
        optimizer::ZoomResult solve = optimizer::solveZoom(rocket, zoom_options);
        min = solve.best_mass;
        index = -1;
        bestDistro = solve.best_units;
        timings.fill_seconds = solve.coarse_timings.fill_seconds;
        timings.tuple_seconds = solve.coarse_timings.tuple_seconds;
        timings.mass_seconds = solve.coarse_timings.mass_seconds + solve.refine_seconds;
        timings.min_seconds = solve.coarse_timings.min_seconds;
        summary.n_combinations = solve.coarse_evaluated + solve.refined_evaluated;
        if (options.print_output) {
            std::cout << "# Combinations (refined): " << solve.refined_evaluated << "\n";
            std::cout << "Refine radius used (fine units): " << solve.radius_used << "\n\n";
        }
    }

    if (!std::isfinite(min)) {
        std::cerr << "No feasible solution found (all distributions invalid).\n";
        setError("No feasible solution found (all distributions invalid).");
        return false;
    }

    auto total_end = Clock::now();
    timings.total_seconds = elapsedSeconds(total_start, total_end);

    summary.timings = timings;
    summary.min_mass = min;
    summary.min_index = index;
    summary.stage_count = static_cast<int>(stage_count);
    summary.precision = effective_precision;
    summary.threads_used = nThreads;
    summary.deltaVmission = rocket.deltaV;
    summary.best_distribution = bestDistro;

    if (options.print_output) {
        std::cout << timings.fill_seconds << " seconds to Fill Array with 1s" << '\n';
        std::cout << timings.tuple_seconds << " seconds to create the tuple of all possible distributions\n\n";
        std::cout << timings.mass_seconds << " seconds to calculate mass of all possible distributions\n\n";
        std::cout << timings.min_seconds << " seconds to find minimum mass" << "\n";
        std::uint64_t min_kg = static_cast<std::uint64_t>(std::ceil(min));
        if (index >= 0) {
            std::cout << min_kg << " kg at Distribution number " << index << "\n";
        }
        else {
            std::cout << min_kg << " kg\n";
        }
        std::cout << (static_cast<double>(min_kg) / 1000.0) << " t" << "\n\n";
        std::cout << "best deltaV distribution is:\n";
        double* newArr = new double[stage_count];
        newArr = ArrToPerc(newArr, bestDistro.data(), static_cast<int>(stage_count), effective_precision);
        outputArr(newArr, static_cast<int>(stage_count));
        std::cout << "\n";
        for (std::size_t i = 0; i < stage_count; i++) {
            std::cout << "Stage " << stage_count - i << ":\t"
                      << ((static_cast<double>(bestDistro[i]) / effective_precision) * rocket.deltaV) << "\n";
        }
        auto warnings = optimizer::buildLowStageDvWarnings(bestDistro.data(),
                                                          stage_count,
                                                          effective_precision,
                                                          rocket.deltaV,
                                                          reporting_min_stage_dv_mps);
        for (const auto& warning : warnings) {
            std::cout << "Warning: " << warning << "\n";
        }
        delete[] newArr;
    }

    if (options.pause_on_exit && options.print_output) {
        std::cout << "\nPress Enter to end programm";
        getchar();
    }

    return true;
}

static BenchmarkTimings sumTimings(const BenchmarkTimings& left, const BenchmarkTimings& right)
{
    BenchmarkTimings out;
    out.fill_seconds = left.fill_seconds + right.fill_seconds;
    out.tuple_seconds = left.tuple_seconds + right.tuple_seconds;
    out.mass_seconds = left.mass_seconds + right.mass_seconds;
    out.min_seconds = left.min_seconds + right.min_seconds;
    out.total_seconds = left.total_seconds + right.total_seconds;
    return out;
}

static BenchmarkTimings divideTimings(const BenchmarkTimings& total, int count)
{
    BenchmarkTimings out;
    out.fill_seconds = total.fill_seconds / count;
    out.tuple_seconds = total.tuple_seconds / count;
    out.mass_seconds = total.mass_seconds / count;
    out.min_seconds = total.min_seconds / count;
    out.total_seconds = total.total_seconds / count;
    return out;
}

static string csvEscape(const string& value)
{
    bool needs_quotes = false;
    for (char c : value) {
        if (c == '"' || c == ',' || c == '\n' || c == '\r') {
            needs_quotes = true;
            break;
        }
    }
    if (!needs_quotes) {
        return value;
    }
    string out;
    out.reserve(value.size() + 2);
    out.push_back('"');
    for (char c : value) {
        if (c == '"') {
            out.push_back('"');
        }
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

static void appendBenchmarkCsv(const string& path,
                               const string& git_head,
                               const string& config_path,
                               const json& config,
                               const OptimizationSummary& summary,
                               const BenchmarkTimings& average_timings,
                               const BenchmarkTimings& worst_timings,
                               const ProgramOptions& options,
                               const string& status)
{
    bool needs_header = true;
    {
        ifstream existing(path);
        if (existing.good()) {
            existing.seekg(0, std::ios::end);
            needs_header = (existing.tellg() == 0);
        }
    }

    ofstream file(path, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Warning: unable to open benchmark CSV at " << path << "\n";
        return;
    }

    if (needs_header) {
        file << "timestamp_utc,benchmark_version,git_head,config_path,rocket_path,engines_path,precision,stages,n_combinations,threads,iterations,warmup,max_total_seconds,fill_avg,tuple_avg,mass_avg,min_avg,total_avg,fill_max,tuple_max,mass_max,min_max,total_max,status\n";
    }

    std::time_t now = std::time(nullptr);
    string rocket_path = config["rocketPath"].get<string>();
    string engines_path = config["enginesPath"].get<string>();

    file << now << ",";
    file << 1 << ",";
    file << csvEscape(git_head) << ",";
    file << csvEscape(config_path) << ",";
    file << csvEscape(rocket_path) << ",";
    file << csvEscape(engines_path) << ",";
    file << summary.precision << ",";
    file << summary.stage_count << ",";
    file << summary.n_combinations << ",";
    file << summary.threads_used << ",";
    file << options.benchmark_iterations << ",";
    file << options.benchmark_warmup << ",";
    if (options.benchmark_max_seconds > 0.0) {
        file << options.benchmark_max_seconds;
    }
    file << ",";
    file << std::fixed << std::setprecision(9)
         << average_timings.fill_seconds << ","
         << average_timings.tuple_seconds << ","
         << average_timings.mass_seconds << ","
         << average_timings.min_seconds << ","
         << average_timings.total_seconds << ","
         << worst_timings.fill_seconds << ","
         << worst_timings.tuple_seconds << ","
         << worst_timings.mass_seconds << ","
         << worst_timings.min_seconds << ","
         << worst_timings.total_seconds << ",";
    file << csvEscape(status);
    file << "\n";
}

static void printDistributionArray(const std::vector<unsigned char>& distro)
{
    std::cout << "[";
    for (std::size_t i = 0; i < distro.size(); i++) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << static_cast<int>(distro[i]);
    }
    std::cout << "]";
}

static void printDistributionDvArray(const std::vector<unsigned char>& distro, int precision, int deltaVmission)
{
    std::cout << "[";
    for (std::size_t i = 0; i < distro.size(); i++) {
        if (i != 0) {
            std::cout << ", ";
        }
        double dv = 0.0;
        if (precision > 0 && deltaVmission > 0) {
            dv = (static_cast<double>(distro[i]) / static_cast<double>(precision)) * static_cast<double>(deltaVmission);
        }
        std::cout << dv;
    }
    std::cout << "]";
}

static int runBenchmark(const json& config, const ProgramOptions& options, const string& config_path)
{
    RunOptions run_options;
    run_options.print_output = false;
    run_options.pause_on_exit = false;
    run_options.force_verbose_off = true;
    run_options.threads_override = options.benchmark_threads;
    run_options.full_grid_mode_override = options.full_grid_mode_override;

    OptimizationSummary summary;
    for (int i = 0; i < options.benchmark_warmup; i++) {
        OptimizationSummary warmup_summary;
        if (!runOptimization(config, run_options, warmup_summary)) {
            return 1;
        }
    }

    BenchmarkTimings total_timings;
    BenchmarkTimings worst_timings;
    bool first_sample = true;

    for (int i = 0; i < options.benchmark_iterations; i++) {
        OptimizationSummary sample_summary;
        if (!runOptimization(config, run_options, sample_summary)) {
            return 1;
        }
        if (first_sample) {
            summary = sample_summary;
            worst_timings = sample_summary.timings;
            first_sample = false;
        }
        total_timings = sumTimings(total_timings, sample_summary.timings);
        worst_timings.total_seconds = std::max(worst_timings.total_seconds, sample_summary.timings.total_seconds);
        worst_timings.fill_seconds = std::max(worst_timings.fill_seconds, sample_summary.timings.fill_seconds);
        worst_timings.tuple_seconds = std::max(worst_timings.tuple_seconds, sample_summary.timings.tuple_seconds);
        worst_timings.mass_seconds = std::max(worst_timings.mass_seconds, sample_summary.timings.mass_seconds);
        worst_timings.min_seconds = std::max(worst_timings.min_seconds, sample_summary.timings.min_seconds);
    }

    BenchmarkTimings average_timings = divideTimings(total_timings, options.benchmark_iterations);
    string git_head = getGitHead();
    bool passed = options.benchmark_max_seconds <= 0.0 || average_timings.total_seconds <= options.benchmark_max_seconds;
    string status = passed ? "pass" : "fail";

    std::cout << std::fixed << std::setprecision(9);
    std::cout << "{\n";
    std::cout << "  \"benchmark_version\": 1,\n";
    std::cout << "  \"git_head\": \"" << jsonEscape(git_head) << "\",\n";
    std::cout << "  \"config_path\": \"" << jsonEscape(config_path) << "\",\n";
    std::cout << "  \"rocket_path\": \"" << jsonEscape(config["rocketPath"].get<string>()) << "\",\n";
    std::cout << "  \"engines_path\": \"" << jsonEscape(config["enginesPath"].get<string>()) << "\",\n";
    std::cout << "  \"precision\": " << summary.precision << ",\n";
    std::cout << "  \"stages\": " << summary.stage_count << ",\n";
    std::cout << "  \"n_combinations\": " << summary.n_combinations << ",\n";
    std::cout << "  \"threads\": " << summary.threads_used << ",\n";
    std::cout << "  \"min_mass\": " << summary.min_mass << ",\n";
    std::cout << "  \"best_distribution_units\": ";
    printDistributionArray(summary.best_distribution);
    std::cout << ",\n";
    std::cout << "  \"best_distribution_dv_mps\": ";
    printDistributionDvArray(summary.best_distribution, summary.precision, summary.deltaVmission);
    std::cout << ",\n";
    std::cout << "  \"iterations\": " << options.benchmark_iterations << ",\n";
    std::cout << "  \"warmup\": " << options.benchmark_warmup << ",\n";
    if (options.benchmark_max_seconds > 0.0) {
        std::cout << "  \"max_total_seconds\": " << options.benchmark_max_seconds << ",\n";
    }
    std::cout << "  \"timings_seconds\": {\n";
    std::cout << "    \"fill_avg\": " << average_timings.fill_seconds << ",\n";
    std::cout << "    \"tuple_avg\": " << average_timings.tuple_seconds << ",\n";
    std::cout << "    \"mass_avg\": " << average_timings.mass_seconds << ",\n";
    std::cout << "    \"min_avg\": " << average_timings.min_seconds << ",\n";
    std::cout << "    \"total_avg\": " << average_timings.total_seconds << ",\n";
    std::cout << "    \"fill_max\": " << worst_timings.fill_seconds << ",\n";
    std::cout << "    \"tuple_max\": " << worst_timings.tuple_seconds << ",\n";
    std::cout << "    \"mass_max\": " << worst_timings.mass_seconds << ",\n";
    std::cout << "    \"min_max\": " << worst_timings.min_seconds << ",\n";
    std::cout << "    \"total_max\": " << worst_timings.total_seconds << "\n";
    std::cout << "  },\n";
    std::cout << "  \"status\": \"" << status << "\"\n";
    std::cout << "}\n";

    appendBenchmarkCsv(options.benchmark_csv_path,
                       git_head,
                       config_path,
                       config,
                       summary,
                       average_timings,
                       worst_timings,
                       options,
                       status);

    return passed ? 0 : 2;
}

struct StrategyBenchmarkResult {
    std::string name;
    std::string status;
    OptimizationSummary summary;
    BenchmarkTimings average_timings;
    BenchmarkTimings worst_timings;
    MemoryCounters peak_memory;
    bool ok = false;
    std::string error;
};

static bool runStrategyBenchmark(const std::string& name,
                                 const json& config,
                                 const ProgramOptions& options,
                                 StrategyBenchmarkResult& out)
{
    out = StrategyBenchmarkResult{};
    out.name = name;

    RunOptions run_options;
    run_options.print_output = false;
    run_options.pause_on_exit = false;
    run_options.force_verbose_off = true;
    run_options.threads_override = options.benchmark_threads;
    run_options.full_grid_mode_override = options.full_grid_mode_override;

    for (int i = 0; i < options.benchmark_warmup; i++) {
        OptimizationSummary warmup_summary;
        std::string error;
        if (!runOptimization(config, run_options, warmup_summary, &error)) {
            out.status = "error";
            out.error = error.empty() ? "Warmup failed" : error;
            return false;
        }
    }

    BenchmarkTimings total_timings;
    BenchmarkTimings worst_timings;
    MemoryCounters worst_memory;
    bool first_sample = true;

    for (int i = 0; i < options.benchmark_iterations; i++) {
        OptimizationSummary sample_summary;
        std::string error;
        MemoryCounters peak;
        bool ok = runWithPeakMemory([&] { return runOptimization(config, run_options, sample_summary, &error); }, peak);
        if (!ok) {
            out.status = "error";
            out.error = error.empty() ? "Optimization failed" : error;
            return false;
        }
        worst_memory.working_set_bytes = std::max(worst_memory.working_set_bytes, peak.working_set_bytes);
        worst_memory.private_bytes = std::max(worst_memory.private_bytes, peak.private_bytes);
        if (first_sample) {
            out.summary = sample_summary;
            worst_timings = sample_summary.timings;
            first_sample = false;
        }
        total_timings = sumTimings(total_timings, sample_summary.timings);
        worst_timings.total_seconds = std::max(worst_timings.total_seconds, sample_summary.timings.total_seconds);
        worst_timings.fill_seconds = std::max(worst_timings.fill_seconds, sample_summary.timings.fill_seconds);
        worst_timings.tuple_seconds = std::max(worst_timings.tuple_seconds, sample_summary.timings.tuple_seconds);
        worst_timings.mass_seconds = std::max(worst_timings.mass_seconds, sample_summary.timings.mass_seconds);
        worst_timings.min_seconds = std::max(worst_timings.min_seconds, sample_summary.timings.min_seconds);
    }

    out.average_timings = divideTimings(total_timings, options.benchmark_iterations);
    out.worst_timings = worst_timings;
    out.peak_memory = worst_memory;

    bool passed = options.benchmark_max_seconds <= 0.0 || out.average_timings.total_seconds <= options.benchmark_max_seconds;
    out.status = passed ? "pass" : "fail";
    out.ok = true;
    return true;
}

static json strategyToJson(const StrategyBenchmarkResult& result)
{
    json out;
    out["name"] = result.name;
    out["status"] = result.status;
    if (!result.error.empty()) {
        out["error"] = result.error;
    }

    out["precision"] = result.summary.precision;
    out["stages"] = result.summary.stage_count;
    out["n_combinations"] = result.summary.n_combinations;
    out["threads"] = result.summary.threads_used;
    out["min_mass"] = result.summary.min_mass;
    out["memory_peak_working_set_bytes"] = result.peak_memory.working_set_bytes;
    out["memory_peak_private_bytes"] = result.peak_memory.private_bytes;

    out["best_distribution_units"] = json::array();
    for (unsigned char v : result.summary.best_distribution) {
        out["best_distribution_units"].push_back(static_cast<int>(v));
    }

    out["best_distribution_dv_mps"] = json::array();
    for (unsigned char v : result.summary.best_distribution) {
        double dv = 0.0;
        if (result.summary.precision > 0 && result.summary.deltaVmission > 0) {
            dv = (static_cast<double>(v) / static_cast<double>(result.summary.precision)) * static_cast<double>(result.summary.deltaVmission);
        }
        out["best_distribution_dv_mps"].push_back(dv);
    }

    out["timings_seconds"] = json::object();
    out["timings_seconds"]["fill_avg"] = result.average_timings.fill_seconds;
    out["timings_seconds"]["tuple_avg"] = result.average_timings.tuple_seconds;
    out["timings_seconds"]["mass_avg"] = result.average_timings.mass_seconds;
    out["timings_seconds"]["min_avg"] = result.average_timings.min_seconds;
    out["timings_seconds"]["total_avg"] = result.average_timings.total_seconds;
    out["timings_seconds"]["fill_max"] = result.worst_timings.fill_seconds;
    out["timings_seconds"]["tuple_max"] = result.worst_timings.tuple_seconds;
    out["timings_seconds"]["mass_max"] = result.worst_timings.mass_seconds;
    out["timings_seconds"]["min_max"] = result.worst_timings.min_seconds;
    out["timings_seconds"]["total_max"] = result.worst_timings.total_seconds;

    return out;
}

static std::string formatUnitsForTable(unsigned char units)
{
    return std::to_string(static_cast<int>(units));
}

static std::string formatDoubleForTable(double value, int precision = 6)
{
    if (!std::isfinite(value)) {
        return "nan";
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

static std::string formatIntForTable(long long value)
{
    return std::to_string(value);
}

static std::string formatUInt64ForTable(std::uint64_t value)
{
    return util::formatThousands(value);
}

static std::string formatMiBForTable(std::uint64_t bytes, int precision = 2)
{
    if (precision < 0) {
        precision = 0;
    }
    double mib = static_cast<double>(bytes) / (1024.0 * 1024.0);
    return util::formatDoubleThousands(mib, precision) + " MiB";
}

static std::string formatMassForTable(double mass_kg)
{
    return util::formatDoubleThousands(mass_kg, 3) + " kg";
}

static void printBenchmarkCompareTable(const std::vector<StrategyBenchmarkResult>& results,
                                      const StrategyBenchmarkResult* baseline)
{
    std::vector<std::string> headers;
    headers.reserve(results.size() + 1);
    headers.push_back("metric");
    for (const auto& r : results) {
        headers.push_back(r.name);
    }

    auto metricRow = [&](const std::string& metric, const std::function<std::string(const StrategyBenchmarkResult&)>& valueFn) {
        std::vector<std::string> row;
        row.reserve(headers.size());
        row.push_back(metric);
        for (const auto& r : results) {
            row.push_back(valueFn(r));
        }
        return row;
    };

    std::vector<std::vector<std::string>> rows;

    rows.push_back(metricRow("status", [](const StrategyBenchmarkResult& r) { return r.status; }));
    rows.push_back(metricRow("error", [](const StrategyBenchmarkResult& r) { return r.error; }));

    rows.push_back(metricRow("precision", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatIntForTable(r.summary.precision);
    }));
    rows.push_back(metricRow("stages", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatIntForTable(r.summary.stage_count);
    }));
    rows.push_back(metricRow("threads", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatIntForTable(r.summary.threads_used);
    }));
    rows.push_back(metricRow("n_combinations", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatUInt64ForTable(static_cast<std::uint64_t>(r.summary.n_combinations));
    }));

    rows.push_back(metricRow("mem_peak_working_set", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error" || r.peak_memory.working_set_bytes == 0) {
            return std::string{};
        }
        return formatMiBForTable(r.peak_memory.working_set_bytes, 2);
    }));
    rows.push_back(metricRow("mem_peak_private", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error" || r.peak_memory.private_bytes == 0) {
            return std::string{};
        }
        return formatMiBForTable(r.peak_memory.private_bytes, 2);
    }));
    rows.push_back(metricRow("min_mass", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatMassForTable(r.summary.min_mass);
    }));

    rows.push_back(metricRow("mass_abs_diff_vs_full_grid", [&](const StrategyBenchmarkResult& r) {
        if (!baseline || r.status == "error" || baseline->status == "error") {
            return std::string{};
        }
        if (r.name == baseline->name) {
            return std::string("0");
        }
        if (!std::isfinite(baseline->summary.min_mass) || !std::isfinite(r.summary.min_mass)) {
            return std::string{};
        }
        return formatDoubleForTable(r.summary.min_mass - baseline->summary.min_mass, 6);
    }));
    rows.push_back(metricRow("mass_rel_diff_vs_full_grid_pct", [&](const StrategyBenchmarkResult& r) {
        if (!baseline || r.status == "error" || baseline->status == "error") {
            return std::string{};
        }
        if (r.name == baseline->name) {
            return std::string("0");
        }
        if (!std::isfinite(baseline->summary.min_mass) || !std::isfinite(r.summary.min_mass) || baseline->summary.min_mass <= 0.0) {
            return std::string{};
        }
        double rel = (r.summary.min_mass / baseline->summary.min_mass) - 1.0;
        return formatDoubleForTable(rel * 100.0, 6);
    }));
    rows.push_back(metricRow("dv_l1_diff_mps_vs_full_grid", [&](const StrategyBenchmarkResult& r) {
        if (!baseline || r.status == "error" || baseline->status == "error") {
            return std::string{};
        }
        if (r.name == baseline->name) {
            return std::string("0");
        }
        if (baseline->summary.stage_count != r.summary.stage_count || baseline->summary.deltaVmission != r.summary.deltaVmission) {
            return std::string{};
        }
        double l1 = 0.0;
        for (std::size_t i = 0; i < baseline->summary.best_distribution.size() && i < r.summary.best_distribution.size(); i++) {
            double dv0 = (baseline->summary.precision > 0) ? (static_cast<double>(baseline->summary.best_distribution[i]) / baseline->summary.precision) * baseline->summary.deltaVmission : 0.0;
            double dv1 = (r.summary.precision > 0) ? (static_cast<double>(r.summary.best_distribution[i]) / r.summary.precision) * r.summary.deltaVmission : 0.0;
            l1 += std::fabs(dv1 - dv0);
        }
        return formatDoubleForTable(l1, 6);
    }));

    rows.push_back(metricRow("total_avg_s", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatDoubleForTable(r.average_timings.total_seconds, 6);
    }));
    rows.push_back(metricRow("total_max_s", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatDoubleForTable(r.worst_timings.total_seconds, 6);
    }));
    rows.push_back(metricRow("min_max_s", [](const StrategyBenchmarkResult& r) {
        if (r.status == "error") {
            return std::string{};
        }
        return formatDoubleForTable(r.worst_timings.min_seconds, 6);
    }));

    int stage_count = 0;
    if (baseline && baseline->status != "error") {
        stage_count = baseline->summary.stage_count;
    }
    else {
        for (const auto& r : results) {
            stage_count = std::max(stage_count, r.summary.stage_count);
        }
    }

    for (int stage_index = 0; stage_index < stage_count; stage_index++) {
        rows.push_back(metricRow("stage" + std::to_string(stage_index + 1) + "_units", [stage_index](const StrategyBenchmarkResult& r) {
            if (r.status == "error") {
                return std::string{};
            }
            if (stage_index < 0 || static_cast<std::size_t>(stage_index) >= r.summary.best_distribution.size()) {
                return std::string{};
            }
            return formatUnitsForTable(r.summary.best_distribution[static_cast<std::size_t>(stage_index)]);
        }));
        rows.push_back(metricRow("stage" + std::to_string(stage_index + 1) + "_dv_mps", [stage_index](const StrategyBenchmarkResult& r) {
            if (r.status == "error") {
                return std::string{};
            }
            if (stage_index < 0 || static_cast<std::size_t>(stage_index) >= r.summary.best_distribution.size()) {
                return std::string{};
            }
            if (r.summary.precision <= 0 || r.summary.deltaVmission <= 0) {
                return std::string{};
            }
            double dv = (static_cast<double>(r.summary.best_distribution[static_cast<std::size_t>(stage_index)]) / r.summary.precision) * r.summary.deltaVmission;
            return formatDoubleForTable(dv, 3);
        }));
    }

    std::cout << util::formatTextTable(headers, rows, 2);
}

static int runBenchmarkCompare(const json& base_config, const ProgramOptions& options, const string& config_path)
{
    ProgramOptions local_options = options;
    local_options.full_grid_mode_override.clear();

    json config = base_config;
    if (!config.contains("verbose")) {
        config["verbose"] = false;
    }
    if (!config.contains("useMultiCore")) {
        config["useMultiCore"] = false;
    }

    int config_precision = 0;
    if (config.contains("precision") && config["precision"].is_number()) {
        config_precision = config["precision"].get<int>();
    }
    int fine_precision = config_precision;
    std::string fine_precision_source = "precision";
    bool has_zoom = config.contains("zoom") && config["zoom"].is_object();
    bool zoom_enabled = false;
    if (has_zoom) {
        const auto& zoom = config["zoom"];
        if (zoom.contains("enabled") && zoom["enabled"].is_boolean()) {
            zoom_enabled = zoom["enabled"].get<bool>();
        }
        if (zoom_enabled && zoom.contains("fine_precision") && zoom["fine_precision"].is_number_integer()) {
            fine_precision = zoom["fine_precision"].get<int>();
            fine_precision_source = "zoom.fine_precision";
        }
    }

    std::string baseline_note;
    if (fine_precision_source == "zoom.fine_precision" && config_precision > 0 && fine_precision > 0 && fine_precision != config_precision) {
        baseline_note = "Note: benchmark-compare baseline 'full_grid' uses precision=" + std::to_string(fine_precision) +
                        " derived from zoom.fine_precision (config precision=" + std::to_string(config_precision) + ").";
    }
    if (!baseline_note.empty()) {
        std::cerr << baseline_note << "\n";
    }

    bool has_budget = config.contains("maxCombinations") || config.contains("targetSeconds");

    std::vector<std::pair<std::string, json>> strategies;

    // Baseline: full grid at fine precision, without budgets.
    {
        json cfg = config;
        cfg.erase("maxCombinations");
        cfg.erase("targetSeconds");
        cfg.erase("zoom");
        if (fine_precision > 0) {
            cfg["precision"] = fine_precision;
        }
        cfg["fullGridMode"] = "materialize";
        strategies.push_back({"full_grid", cfg});

        json streaming_cfg = cfg;
        streaming_cfg["fullGridMode"] = "streaming";
        strategies.push_back({"full_grid_streaming", streaming_cfg});
    }

    // Zoom only (no budgets), if zoom is present in the config.
    if (has_zoom) {
        json cfg = config;
        cfg.erase("maxCombinations");
        cfg.erase("targetSeconds");
        cfg["zoom"]["enabled"] = true;
        if (!cfg["zoom"].contains("fine_precision") || !cfg["zoom"]["fine_precision"].is_number_integer()) {
            cfg["zoom"]["fine_precision"] = fine_precision;
        }
        strategies.push_back({"zoom", cfg});
    }

    // Budget only (no zoom), if budget keys are present.
    if (has_budget) {
        json cfg = config;
        cfg.erase("zoom");
        strategies.push_back({"budget_full", cfg});
    }

    // Budget + zoom, if both are present.
    if (has_budget && has_zoom) {
        json cfg = config;
        cfg["zoom"]["enabled"] = true;
        if (!cfg["zoom"].contains("fine_precision") || !cfg["zoom"]["fine_precision"].is_number_integer()) {
            cfg["zoom"]["fine_precision"] = fine_precision;
        }
        strategies.push_back({"budget_zoom", cfg});
    }

    std::vector<StrategyBenchmarkResult> results;
    results.reserve(strategies.size());

    bool any_error = false;
    bool all_passed = true;

    for (const auto& entry : strategies) {
        StrategyBenchmarkResult r;
        if (!runStrategyBenchmark(entry.first, entry.second, local_options, r)) {
            any_error = true;
            all_passed = false;
        }
        if (r.status != "pass") {
            all_passed = false;
        }
        results.push_back(std::move(r));
    }

    json root;
    root["compare_version"] = 1;
    root["git_head"] = getGitHead();
    root["config_path"] = config_path;
    root["iterations"] = options.benchmark_iterations;
    root["warmup"] = options.benchmark_warmup;
    if (options.benchmark_max_seconds > 0.0) {
        root["max_total_seconds"] = options.benchmark_max_seconds;
    }

    root["strategies"] = json::array();
    for (const auto& r : results) {
        root["strategies"].push_back(strategyToJson(r));
    }

    // Compare against full_grid if available.
    json comparisons = json::array();
    const StrategyBenchmarkResult* baseline = nullptr;
    for (const auto& r : results) {
        if (r.name == "full_grid" && r.status != "error") {
            baseline = &r;
            break;
        }
    }
    if (baseline) {
        for (const auto& r : results) {
            if (r.name == baseline->name || r.status == "error") {
                continue;
            }
            json c;
            c["baseline"] = baseline->name;
            c["name"] = r.name;
            if (std::isfinite(baseline->summary.min_mass) && std::isfinite(r.summary.min_mass) && baseline->summary.min_mass > 0.0) {
                double abs_diff = r.summary.min_mass - baseline->summary.min_mass;
                double rel_diff = (r.summary.min_mass / baseline->summary.min_mass) - 1.0;
                c["mass_abs_diff"] = abs_diff;
                c["mass_rel_diff"] = rel_diff;
            }
            if (baseline->summary.stage_count == r.summary.stage_count && baseline->summary.deltaVmission == r.summary.deltaVmission) {
                double l1 = 0.0;
                for (std::size_t i = 0; i < baseline->summary.best_distribution.size() && i < r.summary.best_distribution.size(); i++) {
                    double dv0 = (baseline->summary.precision > 0) ? (static_cast<double>(baseline->summary.best_distribution[i]) / baseline->summary.precision) * baseline->summary.deltaVmission : 0.0;
                    double dv1 = (r.summary.precision > 0) ? (static_cast<double>(r.summary.best_distribution[i]) / r.summary.precision) * r.summary.deltaVmission : 0.0;
                    l1 += std::fabs(dv1 - dv0);
                }
                c["dv_l1_diff_mps"] = l1;
            }
            comparisons.push_back(std::move(c));
        }
    }
    root["comparisons"] = comparisons;

    std::string comparisons_note;
    if (!baseline) {
        comparisons_note = "Note: baseline strategy 'full_grid' failed; comparisons omitted.";
    }

    if (!baseline_note.empty() || !comparisons_note.empty()) {
        root["notes"] = json::array();
        if (!baseline_note.empty()) {
            root["notes"].push_back(baseline_note);
        }
        if (!comparisons_note.empty()) {
            root["notes"].push_back(comparisons_note);
        }
    }

    root["status"] = all_passed ? "pass" : "fail";

    if (options.benchmark_compare_format == "table") {
        if (!comparisons_note.empty()) {
            std::cout << comparisons_note << "\n";
        }
        if (!comparisons_note.empty()) {
            std::cout << "\n";
        }
        printBenchmarkCompareTable(results, baseline);
    }
    else {
        std::cout << root.dump(2) << "\n";
    }

    if (any_error) {
        return 1;
    }
    return all_passed ? 0 : 2;
}

int main(int argc, char** argv)
{
    ProgramOptions options;
    if (!parseOptions(argc, argv, options)) {
        printUsage();
        return 1;
    }
    if (options.show_help) {
        printUsage();
        return 0;
    }

    string config_path = options.config_path;
    if ((options.benchmark || options.benchmark_compare) && !options.benchmark_config_path.empty()) {
        config_path = options.benchmark_config_path;
    }

    json config;
    try {
        config = readJson(config_path);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    if (options.benchmark_compare) {
        return runBenchmarkCompare(config, options, config_path);
    }

    if (options.benchmark) {
        return runBenchmark(config, options, config_path);
    }

    RunOptions run_options;
    run_options.print_output = true;
    run_options.pause_on_exit = !options.no_prompt;
    run_options.force_verbose_off = false;
    run_options.threads_override = 0;
    run_options.full_grid_mode_override = options.full_grid_mode_override;

    OptimizationSummary summary;
    if (!runOptimization(config, run_options, summary)) {
        return 1;
    }
    return 0;
}
