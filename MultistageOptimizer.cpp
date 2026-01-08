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
#include <chrono>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "rocket_definition/Stage.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "math/math.hpp"
#include "math/massCalc.hpp"
#include "math/validation.hpp"
#include "read_json/read_json.hpp"
#include "util/config_parsing.hpp"

using namespace std;
using MassType = double;

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
};

struct RunOptions {
    bool print_output = true;
    bool pause_on_exit = true;
    bool force_verbose_off = false;
    int threads_override = 0;
};

struct ProgramOptions {
    bool benchmark = false;
    int benchmark_iterations = 1;
    int benchmark_warmup = 0;
    int benchmark_threads = 0;
    double benchmark_max_seconds = -1.0;
    string benchmark_csv_path = "benchmark_results.csv";
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

static void printUsage()
{
    std::cout << "Usage: MultistageOptimizer.exe [options]\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --config <path>                 Config file for normal runs.\n";
    std::cout << "  --no-prompt                      Skip waiting for Enter before exit.\n";
    std::cout << "  --benchmark                      Run performance benchmark.\n";
    std::cout << "  --benchmark-config <path>        Config file for benchmark runs.\n";
    std::cout << "  --benchmark-iterations <count>   Number of benchmark samples.\n";
    std::cout << "  --benchmark-warmup <count>       Warmup runs before sampling.\n";
    std::cout << "  --benchmark-threads <count>      Override thread count.\n";
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

static bool runOptimization(const json& config, const RunOptions& options, OptimizationSummary& summary)
{
    bool verbose = config["verbose"];
    if (options.force_verbose_off) {
        verbose = false;
    }
    bool useMultiCore = config["useMultiCore"];
    unsigned long long maxRAM = 0;
    if (!parseMaxRam(config["maxRAM"], maxRAM)) {
        std::cerr << "Invalid maxRAM value in config.\n";
        return false;
    }
    int precision = config["precision"];

    vector<Engine> engineList;
    Rocket rocket;
    try {
        engineList = importEngines(config.at("enginesPath").get<string>());
        rocket = importRocket(config.at("rocketPath").get<string>(), engineList);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return false;
    }

    const std::size_t stage_count = rocket.stages.size();
    if (stage_count < 1) {
        std::cerr << "Rocket has no stages.\n";
        return false;
    }

    std::uint64_t nCombinations = nCr(precision - 1, static_cast<int>(stage_count) - 1);
    if (nCombinations == 0) {
        std::cerr << "No combinations generated.\n";
        return false;
    }

    unsigned long long distributionsBytes = safeMULT_ULLONG(nCombinations, static_cast<unsigned long long>(stage_count)); // uint8_t
    unsigned long long massesBytes = safeMULT_ULLONG(nCombinations, static_cast<unsigned long long>(sizeof(MassType)));
    unsigned long long requiredBytes = safeAdd_ULLONG(distributionsBytes, massesBytes);
    if (requiredBytes >= maxRAM) {
        std::cerr << "Precision too high for the amount of stages.\n";
        std::cerr << "Maximum precision for current stages and RAM: "
                  << calcMaxPrecUp(precision, static_cast<int>(stage_count), maxRAM, sizeof(MassType)) << "\n";
        return false;
    }

    if (options.print_output) {
        std::cout << "For the number of stages and RAM currently given, the maximum precision can be raised to: "
                  << calcMaxPrecDown(255, static_cast<int>(stage_count), maxRAM, sizeof(MassType)) << "\n\n";
        std::cout << "# Combinations: " << nCombinations << "\n\n";
    }

    int nThreads = resolveThreadCount(useMultiCore, options.threads_override);
    if (options.threads_override > 0) {
        useMultiCore = (nThreads > 1);
    }

    if (nThreads < 1) {
        nThreads = 1;
    }
    if (static_cast<std::uint64_t>(nThreads) > nCombinations) {
        nThreads = static_cast<int>(nCombinations);
        if (nThreads < 1) {
            nThreads = 1;
        }
    }

    unsigned char* distributions = new unsigned char[static_cast<std::size_t>(nCombinations) * stage_count];
    MassType* massDistributions = new MassType[static_cast<std::size_t>(nCombinations)];

    auto total_start = Clock::now();

    auto fill_start = Clock::now();
    std::fill_n(distributions, nCombinations * rocket.stages.size(), 1);
    auto fill_end = Clock::now();

    auto tuple_start = Clock::now();
    std::uint64_t generated = createTuple(distributions, rocket.stages.size(), nCombinations, precision);
    if (generated != nCombinations) {
        std::cerr << "Tuple generation failed (generated " << generated << " of " << nCombinations << ").\n";
        delete[] distributions;
        delete[] massDistributions;
        return false;
    }
    auto tuple_end = Clock::now();

    if (verbose) {
        if (!validateDistributions(distributions, stage_count, nCombinations, precision, 100000)) {
            std::cerr << "Distribution validation failed.\n";
            delete[] distributions;
            delete[] massDistributions;
            return false;
        }
    }

    std::fill_n(massDistributions, nCombinations, 0.0);

    auto mass_start = Clock::now();
    if (options.print_output) {
        if (useMultiCore) {
            cout << "using " << nThreads << " threads\n";
        }
        else {
            cout << "using 1 thread\n";
        }
    }

    if (nThreads == 1) {
        distMass(verbose, 0, nCombinations, nCombinations, rocket, massDistributions, distributions, precision);
    }
    else {
        std::uint64_t threadBlock = nCombinations / static_cast<std::uint64_t>(nThreads);
        std::vector<std::thread> threads;
        threads.reserve(static_cast<std::size_t>(nThreads));

        for (int t = 0; t < nThreads; t++) {
            std::uint64_t begin = static_cast<std::uint64_t>(t) * threadBlock;
            std::uint64_t end = (t == nThreads - 1) ? nCombinations : begin + threadBlock;
            threads.push_back(std::thread([verbose, begin, end, nCombinations, &rocket, massDistributions, distributions, precision] {
                distMass(verbose, begin, end, nCombinations, rocket, massDistributions, distributions, precision);
            }));
        }

        for (auto& th : threads) {
            th.join();
        }
    }
    auto mass_end = Clock::now();

    auto min_start = Clock::now();
    MassType min = std::numeric_limits<MassType>::infinity();
    std::uint64_t index = 0;
    unsigned char* bestDistro = new unsigned char[stage_count];
    for (std::uint64_t i = 0; i < nCombinations; i++) {
        if (massDistributions[static_cast<std::size_t>(i)] < min) {
            min = massDistributions[static_cast<std::size_t>(i)];
            index = i;
            for (std::size_t j = 0; j < stage_count; j++) {
                bestDistro[j] = distributions[j * nCombinations + i];
            }
        }
    }
    if (!std::isfinite(min)) {
        std::cerr << "No feasible solution found (all distributions invalid).\n";
        delete[] distributions;
        delete[] massDistributions;
        delete[] bestDistro;
        return false;
    }
    auto min_end = Clock::now();
    auto total_end = Clock::now();

    BenchmarkTimings timings;
    timings.fill_seconds = elapsedSeconds(fill_start, fill_end);
    timings.tuple_seconds = elapsedSeconds(tuple_start, tuple_end);
    timings.mass_seconds = elapsedSeconds(mass_start, mass_end);
    timings.min_seconds = elapsedSeconds(min_start, min_end);
    timings.total_seconds = elapsedSeconds(total_start, total_end);

    summary.timings = timings;
    summary.n_combinations = nCombinations;
    summary.min_mass = min;
    summary.min_index = static_cast<long long>(index);
    summary.stage_count = static_cast<int>(stage_count);
    summary.precision = precision;
    summary.threads_used = nThreads;

    if (options.print_output) {
        std::cout << timings.fill_seconds << " seconds to Fill Array with 1s" << '\n';
        std::cout << timings.tuple_seconds << " seconds to create the tuple of all possible distributions\n\n";
        std::cout << timings.mass_seconds << " seconds to calculate mass of all possible distributions\n\n";
        std::cout << timings.min_seconds << " seconds to find minimum mass" << "\n";
        std::uint64_t min_kg = static_cast<std::uint64_t>(std::ceil(min));
        std::cout << min_kg << " kg at Distribution number " << index << "\n";
        std::cout << (static_cast<double>(min_kg) / 1000.0) << " t" << "\n\n";
        std::cout << "best deltaV distribution is:\n";
        double* newArr = new double[stage_count];
        newArr = ArrToPerc(newArr, bestDistro, static_cast<int>(stage_count), precision);
        outputArr(newArr, static_cast<int>(stage_count));
        std::cout << "\n";
        for (std::size_t i = 0; i < stage_count; i++) {
            std::cout << "Stage " << stage_count - i << ":\t"
                      << ((static_cast<double>(bestDistro[i]) / precision) * rocket.deltaV) << "\n";
        }
        delete[] newArr;
    }

    delete[] distributions;
    delete[] massDistributions;
    delete[] bestDistro;

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

static int runBenchmark(const json& config, const ProgramOptions& options, const string& config_path)
{
    RunOptions run_options;
    run_options.print_output = false;
    run_options.pause_on_exit = false;
    run_options.force_verbose_off = true;
    run_options.threads_override = options.benchmark_threads;

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
    if (options.benchmark && !options.benchmark_config_path.empty()) {
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

    if (options.benchmark) {
        return runBenchmark(config, options, config_path);
    }

    RunOptions run_options;
    run_options.print_output = true;
    run_options.pause_on_exit = !options.no_prompt;
    run_options.force_verbose_off = false;
    run_options.threads_override = 0;

    OptimizationSummary summary;
    if (!runOptimization(config, run_options, summary)) {
        return 1;
    }
    return 0;
}
