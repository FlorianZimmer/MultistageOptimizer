#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "math/massCalc.hpp"
#include "math/math.hpp"
#include "math/validation.hpp"
#include "optimizer/bounded_compositions.hpp"
#include "rocket_definition/Rocket.hpp"

namespace optimizer {

struct SolveTimings {
    double fill_seconds = 0.0;
    double tuple_seconds = 0.0;
    double mass_seconds = 0.0;
    double min_seconds = 0.0;
};

struct Candidate {
    double mass = std::numeric_limits<double>::infinity();
    std::uint64_t index = 0;
    std::vector<unsigned char> units;
};

struct FullGridResult {
    std::uint64_t combinations = 0;
    double best_mass = std::numeric_limits<double>::infinity();
    std::uint64_t best_index = 0;
    std::vector<unsigned char> best_units;
    std::vector<Candidate> best_k;
    SolveTimings timings;
};

enum class FullGridMode {
    Materialize,
    Streaming,
};

struct RocketEvalTables {
    int precision = 0;
    int deltaV = 0;
    double g_body = 0.0;
    std::vector<double> inert_mass_rest;
    std::vector<double> inert_mass_fract;
    std::vector<double> one_minus_inert_fract;
    std::vector<double> inert_coeff;
    std::vector<std::vector<double>> ratio_tables;
    bool valid = false;
};

inline RocketEvalTables buildRocketEvalTables(const Rocket& rocket, int precision)
{
    RocketEvalTables out;
    out.precision = precision;
    out.deltaV = rocket.deltaV;
    out.g_body = rocket.g_body;

    const std::size_t stages = rocket.stages.size();
    if (stages == 0) {
        return out;
    }
    if (precision <= 0 || precision > 255) {
        return out;
    }
    if (!std::isfinite(out.g_body) || out.g_body <= 0.0) {
        return out;
    }
    if (out.deltaV < 0) {
        return out;
    }

    out.inert_mass_rest.resize(stages);
    out.inert_mass_fract.resize(stages);
    out.one_minus_inert_fract.resize(stages);
    out.inert_coeff.resize(stages);
    out.ratio_tables.resize(stages);

    for (std::size_t j = 0; j < stages; j++) {
        out.inert_mass_rest[j] = rocket.stages[j].inert_mass_rest;
        out.inert_mass_fract[j] = rocket.stages[j].inert_mass_fract;
        double isp_vac = rocket.stages[j].engine.isp_vac;
        if (!std::isfinite(out.inert_mass_rest[j]) || out.inert_mass_rest[j] < 0.0) {
            return out;
        }
        if (!std::isfinite(out.inert_mass_fract[j]) || out.inert_mass_fract[j] <= 0.0 || out.inert_mass_fract[j] >= 1.0) {
            return out;
        }
        if (!std::isfinite(isp_vac) || isp_vac <= 0.0) {
            return out;
        }

        out.one_minus_inert_fract[j] = 1.0 - out.inert_mass_fract[j];
        out.inert_coeff[j] = out.inert_mass_fract[j] / out.one_minus_inert_fract[j];

        out.ratio_tables[j].resize(static_cast<std::size_t>(precision) + 1);
        double coeff = static_cast<double>(out.deltaV) / (static_cast<double>(precision) * out.g_body * isp_vac);
        for (int v = 0; v <= precision; v++) {
            out.ratio_tables[j][static_cast<std::size_t>(v)] = std::exp(static_cast<double>(v) * coeff);
        }
    }

    out.valid = true;
    return out;
}

inline double evaluateDistributionMass(const Rocket& rocket, const RocketEvalTables& tables, const unsigned char* dv_units)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    if (!tables.valid || dv_units == nullptr) {
        return inf;
    }

    const std::size_t stages = rocket.stages.size();
    double current = inf;
    for (std::size_t j = 0; j < stages; j++) {
        unsigned int units = dv_units[j];
        if (units < 1 || units > static_cast<unsigned int>(tables.precision)) {
            return inf;
        }
        double ratio = tables.ratio_tables[j][units];
        if (j == 0) {
            current = calcMassFromRatio(rocket.stages[j].payload,
                                        tables.inert_mass_fract[j],
                                        tables.one_minus_inert_fract[j],
                                        tables.inert_coeff[j],
                                        ratio,
                                        tables.inert_mass_rest[j]);
        }
        else {
            current = calcMassFromRatio(current,
                                        tables.inert_mass_fract[j],
                                        tables.one_minus_inert_fract[j],
                                        tables.inert_coeff[j],
                                        ratio,
                                        tables.inert_mass_rest[j]);
        }
        if (!std::isfinite(current)) {
            return inf;
        }
    }
    return current;
}

inline double evaluateDistributionMass(const Rocket& rocket, const RocketEvalTables& tables, const int* dv_units)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    if (!tables.valid || dv_units == nullptr) {
        return inf;
    }

    const std::size_t stages = rocket.stages.size();
    double current = inf;
    for (std::size_t j = 0; j < stages; j++) {
        int units = dv_units[j];
        if (units < 1 || units > tables.precision) {
            return inf;
        }
        double ratio = tables.ratio_tables[j][static_cast<std::size_t>(units)];
        if (j == 0) {
            current = calcMassFromRatio(rocket.stages[j].payload,
                                        tables.inert_mass_fract[j],
                                        tables.one_minus_inert_fract[j],
                                        tables.inert_coeff[j],
                                        ratio,
                                        tables.inert_mass_rest[j]);
        }
        else {
            current = calcMassFromRatio(current,
                                        tables.inert_mass_fract[j],
                                        tables.one_minus_inert_fract[j],
                                        tables.inert_coeff[j],
                                        ratio,
                                        tables.inert_mass_rest[j]);
        }
        if (!std::isfinite(current)) {
            return inf;
        }
    }
    return current;
}

inline bool betterCandidate(double lhs_mass,
                            const std::vector<unsigned char>& lhs_units,
                            double rhs_mass,
                            const std::vector<unsigned char>& rhs_units,
                            double rel_eps = 1e-12)
{
    if (lhs_mass < rhs_mass) {
        return true;
    }
    if (!std::isfinite(lhs_mass) || !std::isfinite(rhs_mass)) {
        return false;
    }
    double diff = std::fabs(lhs_mass - rhs_mass);
    double scale = std::max(1.0, std::max(std::fabs(lhs_mass), std::fabs(rhs_mass)));
    if (diff > rel_eps * scale) {
        return false;
    }
    return lhs_units < rhs_units;
}

inline FullGridResult solveFullGridStreaming(const Rocket& rocket,
                                            int precision,
                                            bool verbose,
                                            int threads,
                                            std::size_t keep_top_k = 1)
{
    FullGridResult result;

    using Clock = std::chrono::steady_clock;
    auto elapsedSeconds = [](Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double>(end - start).count();
    };

    const std::size_t stages = rocket.stages.size();
    if (stages == 0) {
        return result;
    }
    if (precision <= 0 || precision > 255) {
        return result;
    }

    const std::uint64_t nCombinations = nCr(precision - 1, static_cast<int>(stages) - 1);
    result.combinations = nCombinations;
    if (nCombinations == 0) {
        return result;
    }

    if (threads < 1) {
        threads = 1;
    }
    if (static_cast<std::uint64_t>(threads) > nCombinations) {
        threads = static_cast<int>(nCombinations);
        if (threads < 1) {
            threads = 1;
        }
    }

    auto mass_start = Clock::now();
    RocketEvalTables tables = buildRocketEvalTables(rocket, precision);
    if (!tables.valid) {
        return result;
    }

    if (keep_top_k < 1) {
        keep_top_k = 1;
    }

    auto lessByMassIndex = [](const Candidate& left, const Candidate& right) {
        if (left.mass != right.mass) {
            return left.mass < right.mass;
        }
        return left.index < right.index;
    };
    auto betterByMassIndex = [](double mass, std::uint64_t index, const Candidate& other) {
        if (mass != other.mass) {
            return mass < other.mass;
        }
        return index < other.index;
    };

    if (stages == 1) {
        Candidate c;
        c.index = 0;
        c.units.assign(1, static_cast<unsigned char>(precision));
        c.mass = calcMassFromRatio(rocket.stages[0].payload,
                                   tables.inert_mass_fract[0],
                                   tables.one_minus_inert_fract[0],
                                   tables.inert_coeff[0],
                                   tables.ratio_tables[0][static_cast<std::size_t>(precision)],
                                   tables.inert_mass_rest[0]);
        if (!std::isfinite(c.mass)) {
            return result;
        }
        result.best_k.clear();
        result.best_k.push_back(c);
        result.best_mass = c.mass;
        result.best_index = c.index;
        result.best_units = c.units;
        result.timings.fill_seconds = 0.0;
        result.timings.tuple_seconds = 0.0;
        result.timings.mass_seconds = 0.0;
        result.timings.min_seconds = 0.0;
        return result;
    }

    const int max_first = precision - static_cast<int>(stages) + 1;
    if (max_first < 1) {
        return result;
    }

    std::vector<std::uint64_t> prefix(static_cast<std::size_t>(max_first) + 2, 0);
    std::vector<std::uint64_t> count_per_first(static_cast<std::size_t>(max_first) + 1, 0);
    std::uint64_t running = 0;
    const int rest_stages = static_cast<int>(stages) - 1;
    for (int u0 = 1; u0 <= max_first; u0++) {
        prefix[static_cast<std::size_t>(u0)] = running;
        int rest_sum = precision - u0;
        std::uint64_t count = 0;
        if (rest_stages == 1) {
            count = 1;
        }
        else {
            unsigned long long count_ull = nCr(rest_sum - 1, rest_stages - 1);
            count = (count_ull == ULLONG_MAX) ? std::numeric_limits<std::uint64_t>::max() : static_cast<std::uint64_t>(count_ull);
        }
        count_per_first[static_cast<std::size_t>(u0)] = count;
        running = safeAdd_ULLONG(running, count);
    }

    if (threads > max_first) {
        threads = max_first;
    }
    if (threads < 1) {
        threads = 1;
    }

    struct FirstRange {
        int begin = 1;
        int end = 1; // exclusive
    };

    std::vector<FirstRange> ranges;
    ranges.reserve(static_cast<std::size_t>(threads));

    std::uint64_t target = nCombinations / static_cast<std::uint64_t>(threads);
    if (target == 0) {
        target = 1;
    }

    int current_begin = 1;
    std::uint64_t bucket = 0;
    for (int u0 = 1; u0 <= max_first; u0++) {
        std::uint64_t c = count_per_first[static_cast<std::size_t>(u0)];
        if (ranges.size() + 1 < static_cast<std::size_t>(threads) && bucket > 0 && bucket + c > target) {
            ranges.push_back(FirstRange{current_begin, u0});
            current_begin = u0;
            bucket = 0;
        }
        bucket = safeAdd_ULLONG(bucket, c);
    }
    ranges.push_back(FirstRange{current_begin, max_first + 1});

    std::vector<std::vector<Candidate>> top_per_thread(ranges.size());

    auto maybe_insert = [&](std::vector<Candidate>& top, double mass, std::uint64_t index, const unsigned char* units) {
        if (!std::isfinite(mass)) {
            return;
        }
        if (top.size() >= keep_top_k && !betterByMassIndex(mass, index, top.back())) {
            return;
        }
        Candidate candidate;
        candidate.mass = mass;
        candidate.index = index;
        candidate.units.assign(stages, 0);
        for (std::size_t i = 0; i < stages; i++) {
            candidate.units[i] = units[i];
        }
        auto pos = std::lower_bound(top.begin(), top.end(), candidate, lessByMassIndex);
        top.insert(pos, std::move(candidate));
        if (top.size() > keep_top_k) {
            top.pop_back();
        }
    };

    auto worker = [&](std::size_t slot, int u0_begin, int u0_end, bool do_verbose) {
        constexpr double inf = std::numeric_limits<double>::infinity();
        auto& local_top = top_per_thread[slot];
        local_top.clear();
        local_top.reserve(keep_top_k);

        std::vector<unsigned char> units_u8;
        if (do_verbose) {
            units_u8.assign(stages, 0);
        }
        double verbose_min = inf;

        std::vector<int> minUnits(static_cast<std::size_t>(rest_stages), 1);
        std::vector<int> maxUnits(static_cast<std::size_t>(rest_stages), 1);

        for (int u0 = u0_begin; u0 < u0_end; u0++) {
            int rest_sum = precision - u0;
            std::uint64_t index_base = prefix[static_cast<std::size_t>(u0)];

            unsigned char u0_u8 = static_cast<unsigned char>(u0);
            double current0 = calcMassFromRatio(rocket.stages[0].payload,
                                                tables.inert_mass_fract[0],
                                                tables.one_minus_inert_fract[0],
                                                tables.inert_coeff[0],
                                                tables.ratio_tables[0][u0_u8],
                                                tables.inert_mass_rest[0]);
            if (!std::isfinite(current0)) {
                continue;
            }

            if (rest_stages == 1) {
                int u1 = rest_sum;
                unsigned char u1_u8 = static_cast<unsigned char>(u1);
                double current = calcMassFromRatio(current0,
                                                  tables.inert_mass_fract[1],
                                                  tables.one_minus_inert_fract[1],
                                                  tables.inert_coeff[1],
                                                  tables.ratio_tables[1][u1_u8],
                                                  tables.inert_mass_rest[1]);
                if (!std::isfinite(current)) {
                    continue;
                }

                if (do_verbose) {
                    units_u8[0] = u0_u8;
                    units_u8[1] = u1_u8;
                    if (current < verbose_min) {
                        verbose_min = current;
                        std::cout << "Distribution number: " << index_base << "\t" << current << "\n";
                        for (std::size_t j = 0; j < stages; j++) {
                            std::cout << "Stage " << stages - j << ": " << static_cast<int>(units_u8[j]) << "\t";
                        }
                        std::cout << "\n\n";
                    }
                }

                unsigned char scratch[2] = {u0_u8, u1_u8};
                maybe_insert(local_top, current, index_base, scratch);
                continue;
            }

            int global_max = rest_sum - (rest_stages - 1);
            maxUnits.assign(static_cast<std::size_t>(rest_stages), global_max);

            std::uint64_t local = 0;
            enumerateBoundedCompositions(minUnits, maxUnits, rest_sum, [&](const std::vector<int>& rest) {
                std::uint64_t index = index_base + local;
                local++;

                double current = current0;
                for (int s = 1; s < static_cast<int>(stages); s++) {
                    int u = rest[static_cast<std::size_t>(s - 1)];
                    unsigned char u_u8 = static_cast<unsigned char>(u);
                    double ratio = tables.ratio_tables[static_cast<std::size_t>(s)][u_u8];
                    current = calcMassFromRatio(current,
                                                tables.inert_mass_fract[static_cast<std::size_t>(s)],
                                                tables.one_minus_inert_fract[static_cast<std::size_t>(s)],
                                                tables.inert_coeff[static_cast<std::size_t>(s)],
                                                ratio,
                                                tables.inert_mass_rest[static_cast<std::size_t>(s)]);
                    if (!std::isfinite(current)) {
                        return;
                    }
                    if (do_verbose) {
                        units_u8[static_cast<std::size_t>(s)] = u_u8;
                    }
                }

                if (do_verbose) {
                    units_u8[0] = u0_u8;
                    if (current < verbose_min) {
                        verbose_min = current;
                        std::cout << "Distribution number: " << index << "\t" << current << "\n";
                        for (std::size_t j = 0; j < stages; j++) {
                            std::cout << "Stage " << stages - j << ": " << static_cast<int>(units_u8[j]) << "\t";
                        }
                        std::cout << "\n\n";
                    }
                }

                if (local_top.size() >= keep_top_k && !betterByMassIndex(current, index, local_top.back())) {
                    return;
                }

                std::vector<unsigned char> scratch_units;
                scratch_units.assign(stages, 0);
                scratch_units[0] = u0_u8;
                for (int s = 1; s < static_cast<int>(stages); s++) {
                    scratch_units[static_cast<std::size_t>(s)] = static_cast<unsigned char>(rest[static_cast<std::size_t>(s - 1)]);
                }
                maybe_insert(local_top, current, index, scratch_units.data());
            });
        }
    };

    if (ranges.size() == 1) {
        worker(0, ranges[0].begin, ranges[0].end, verbose);
    }
    else {
        std::vector<std::thread> pool;
        pool.reserve(ranges.size());
        for (std::size_t t = 0; t < ranges.size(); t++) {
            pool.push_back(std::thread([&, t] {
                worker(t, ranges[t].begin, ranges[t].end, false);
            }));
        }
        for (auto& th : pool) {
            th.join();
        }
    }
    auto mass_end = Clock::now();

    auto min_start = Clock::now();
    std::vector<Candidate> top;
    top.reserve(std::min<std::size_t>(keep_top_k * top_per_thread.size(), static_cast<std::size_t>(std::numeric_limits<int>::max())));
    for (const auto& local : top_per_thread) {
        for (const auto& entry : local) {
            top.push_back(entry);
        }
    }
    std::sort(top.begin(), top.end(), lessByMassIndex);
    if (top.size() > keep_top_k) {
        top.resize(keep_top_k);
    }

    result.best_k = top;

    if (!result.best_k.empty()) {
        result.best_mass = result.best_k.front().mass;
        result.best_index = result.best_k.front().index;
        result.best_units = result.best_k.front().units;
    }

    auto min_end = Clock::now();

    result.timings.fill_seconds = 0.0;
    result.timings.tuple_seconds = 0.0;
    result.timings.mass_seconds = elapsedSeconds(mass_start, mass_end);
    result.timings.min_seconds = elapsedSeconds(min_start, min_end);

    return result;
}

inline FullGridResult solveFullGridMaterialized(const Rocket& rocket,
                                               int precision,
                                               bool verbose,
                                               int threads,
                                               std::size_t keep_top_k = 1)
{
    FullGridResult result;

    using Clock = std::chrono::steady_clock;
    auto elapsedSeconds = [](Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double>(end - start).count();
    };

    const std::size_t stages = rocket.stages.size();
    if (stages == 0) {
        return result;
    }
    if (precision <= 0 || precision > 255) {
        return result;
    }

    const std::uint64_t nCombinations = nCr(precision - 1, static_cast<int>(stages) - 1);
    result.combinations = nCombinations;
    if (nCombinations == 0) {
        return result;
    }

    if (threads < 1) {
        threads = 1;
    }
    if (static_cast<std::uint64_t>(threads) > nCombinations) {
        threads = static_cast<int>(nCombinations);
        if (threads < 1) {
            threads = 1;
        }
    }

    if (keep_top_k < 1) {
        keep_top_k = 1;
    }

    std::vector<unsigned char> distributions(static_cast<std::size_t>(nCombinations) * stages, 0);

    auto fill_start = Clock::now();
    std::fill(distributions.begin(), distributions.end(), static_cast<unsigned char>(1));
    auto fill_end = Clock::now();

    auto tuple_start = Clock::now();
    std::uint64_t generated = createTuple(distributions.data(), stages, nCombinations, precision);
    auto tuple_end = Clock::now();
    if (generated != nCombinations) {
        return result;
    }

    if (verbose) {
        if (!validateDistributions(distributions.data(), stages, nCombinations, precision, 100000)) {
            return result;
        }
    }

    auto mass_start = Clock::now();
    RocketEvalTables tables = buildRocketEvalTables(rocket, precision);
    if (!tables.valid) {
        return result;
    }

    struct IndexCandidate {
        double mass = std::numeric_limits<double>::infinity();
        std::uint64_t index = 0;
    };
    auto lessByMassIndex = [](const IndexCandidate& left, const IndexCandidate& right) {
        if (left.mass != right.mass) {
            return left.mass < right.mass;
        }
        return left.index < right.index;
    };

    std::vector<std::vector<IndexCandidate>> top_per_thread(static_cast<std::size_t>(threads));

    auto scan_range = [&](std::uint64_t begin, std::uint64_t end, std::vector<IndexCandidate>& top, bool do_verbose) {
        constexpr double inf = std::numeric_limits<double>::infinity();
        top.clear();
        top.reserve(keep_top_k);

        std::vector<unsigned char> dv_units;
        if (do_verbose) {
            dv_units.assign(stages, 0);
        }
        double verbose_min = inf;

        for (std::uint64_t i = begin; i < end; i++) {
            double current = inf;
            for (std::size_t j = 0; j < stages; j++) {
                unsigned char units = distributions[j * nCombinations + i];
                if (units < 1 || units > static_cast<unsigned char>(precision)) {
                    current = inf;
                    break;
                }
                if (do_verbose) {
                    dv_units[j] = units;
                }
                double ratio = tables.ratio_tables[j][units];
                if (j == 0) {
                    current = calcMassFromRatio(rocket.stages[0].payload,
                                                tables.inert_mass_fract[0],
                                                tables.one_minus_inert_fract[0],
                                                tables.inert_coeff[0],
                                                ratio,
                                                tables.inert_mass_rest[0]);
                }
                else {
                    current = calcMassFromRatio(current,
                                                tables.inert_mass_fract[j],
                                                tables.one_minus_inert_fract[j],
                                                tables.inert_coeff[j],
                                                ratio,
                                                tables.inert_mass_rest[j]);
                }
                if (!std::isfinite(current)) {
                    current = inf;
                    break;
                }
            }
            if (!std::isfinite(current)) {
                continue;
            }

            if (do_verbose && current < verbose_min) {
                verbose_min = current;
                std::cout << "Distribution number: " << i << "\t" << current << "\n";
                for (std::size_t j = 0; j < stages; j++) {
                    std::cout << "Stage " << stages - j << ": " << static_cast<int>(dv_units[j]) << "\t";
                }
                std::cout << "\n\n";
            }

            IndexCandidate candidate;
            candidate.mass = current;
            candidate.index = i;

            if (top.size() < keep_top_k) {
                auto pos = std::lower_bound(top.begin(), top.end(), candidate, lessByMassIndex);
                top.insert(pos, candidate);
                continue;
            }

            if (top.empty()) {
                top.push_back(candidate);
                continue;
            }

            if (lessByMassIndex(candidate, top.back())) {
                auto pos = std::lower_bound(top.begin(), top.end(), candidate, lessByMassIndex);
                top.insert(pos, candidate);
                if (top.size() > keep_top_k) {
                    top.pop_back();
                }
            }
        }
    };

    if (threads == 1) {
        scan_range(0, nCombinations, top_per_thread[0], verbose);
    }
    else {
        std::uint64_t threadBlock = nCombinations / static_cast<std::uint64_t>(threads);
        std::vector<std::thread> pool;
        pool.reserve(static_cast<std::size_t>(threads));
        for (int t = 0; t < threads; t++) {
            std::uint64_t begin = static_cast<std::uint64_t>(t) * threadBlock;
            std::uint64_t end = (t == threads - 1) ? nCombinations : begin + threadBlock;
            pool.push_back(std::thread([&, t, begin, end] {
                scan_range(begin, end, top_per_thread[static_cast<std::size_t>(t)], false);
            }));
        }
        for (auto& th : pool) {
            th.join();
        }
    }
    auto mass_end = Clock::now();

    auto min_start = Clock::now();
    std::vector<IndexCandidate> top;
    top.reserve(std::min<std::size_t>(keep_top_k * static_cast<std::size_t>(threads), static_cast<std::size_t>(std::numeric_limits<int>::max())));
    for (const auto& local : top_per_thread) {
        for (const auto& entry : local) {
            top.push_back(entry);
        }
    }
    std::sort(top.begin(), top.end(), lessByMassIndex);
    if (top.size() > keep_top_k) {
        top.resize(keep_top_k);
    }

    result.best_k.clear();
    result.best_k.reserve(top.size());
    for (const auto& entry : top) {
        Candidate c;
        c.mass = entry.mass;
        c.index = entry.index;
        c.units.assign(stages, 0);
        for (std::size_t j = 0; j < stages; j++) {
            c.units[j] = distributions[j * nCombinations + entry.index];
        }
        result.best_k.push_back(std::move(c));
    }

    if (!result.best_k.empty()) {
        result.best_mass = result.best_k.front().mass;
        result.best_index = result.best_k.front().index;
        result.best_units = result.best_k.front().units;
    }

    auto min_end = Clock::now();

    result.timings.fill_seconds = elapsedSeconds(fill_start, fill_end);
    result.timings.tuple_seconds = elapsedSeconds(tuple_start, tuple_end);
    result.timings.mass_seconds = elapsedSeconds(mass_start, mass_end);
    result.timings.min_seconds = elapsedSeconds(min_start, min_end);

    return result;
}

inline FullGridResult solveFullGrid(const Rocket& rocket,
                                   int precision,
                                   bool verbose,
                                   int threads,
                                   std::size_t keep_top_k = 1,
                                   FullGridMode mode = FullGridMode::Materialize)
{
    if (mode == FullGridMode::Streaming) {
        return solveFullGridStreaming(rocket, precision, verbose, threads, keep_top_k);
    }
    return solveFullGridMaterialized(rocket, precision, verbose, threads, keep_top_k);
}

struct RefineResult {
    double best_mass = std::numeric_limits<double>::infinity();
    std::vector<unsigned char> best_units;
    std::uint64_t evaluated = 0;
    int used_radius = 0;
};

inline RefineResult refineAroundCenter(const Rocket& rocket,
                                      int finePrecision,
                                      const std::vector<int>& centerUnits,
                                      int initialRadius,
                                      std::uint64_t maxEvaluations = std::numeric_limits<std::uint64_t>::max())
{
    RefineResult out;
    const std::size_t stages = rocket.stages.size();
    if (stages == 0 || centerUnits.size() != stages) {
        return out;
    }
    if (finePrecision <= 0 || finePrecision > 255) {
        return out;
    }
    if (initialRadius < 0) {
        initialRadius = 0;
    }

    RocketEvalTables tables = buildRocketEvalTables(rocket, finePrecision);
    if (!tables.valid) {
        return out;
    }

    std::vector<int> minUnits;
    std::vector<int> maxUnits;

    int radius = initialRadius;
    std::uint64_t count = 0;
    const std::uint64_t count_cap = (maxEvaluations == std::numeric_limits<std::uint64_t>::max())
                                        ? std::numeric_limits<std::uint64_t>::max()
                                        : (maxEvaluations + 1);
    for (;;) {
        computeBoundsAroundCenter(centerUnits, finePrecision, radius, minUnits, maxUnits);
        count = countBoundedCompositions(minUnits, maxUnits, finePrecision, count_cap);
        if (count <= maxEvaluations) {
            break;
        }
        if (radius == 0) {
            count = 1;
            break;
        }
        radius--;
    }
    out.used_radius = radius;

    out.best_units.assign(stages, 0);
    std::vector<int> current(stages, 0);

    auto consider = [&](const std::vector<int>& units) {
        for (std::size_t i = 0; i < stages; i++) {
            current[i] = units[i];
        }
        double mass = evaluateDistributionMass(rocket, tables, current.data());
        if (!std::isfinite(mass)) {
            return;
        }
        std::vector<unsigned char> u8(stages, 0);
        for (std::size_t i = 0; i < stages; i++) {
            u8[i] = static_cast<unsigned char>(units[i]);
        }
        if (betterCandidate(mass, u8, out.best_mass, out.best_units)) {
            out.best_mass = mass;
            out.best_units = std::move(u8);
        }
    };

    enumerateBoundedCompositions(minUnits, maxUnits, finePrecision, consider, maxEvaluations);
    out.evaluated = countBoundedCompositions(minUnits, maxUnits, finePrecision);

    return out;
}

struct ZoomOptions {
    int coarse_precision = 0;
    int fine_precision = 0;
    int window_coarse_steps = 0;
    std::size_t top_k = 1;
    int threads = 1;
    FullGridMode full_grid_mode = FullGridMode::Materialize;
    bool verbose = false;
    std::uint64_t max_evaluations = std::numeric_limits<std::uint64_t>::max();
};

struct ZoomResult {
    int coarse_precision_used = 0;
    int fine_precision_used = 0;
    int radius_used = 0;
    std::uint64_t coarse_evaluated = 0;
    std::uint64_t refined_evaluated = 0;
    double best_mass = std::numeric_limits<double>::infinity();
    std::vector<unsigned char> best_units;
    std::vector<unsigned char> coarse_best_units;
    SolveTimings coarse_timings;
    double refine_seconds = 0.0;
};

inline ZoomResult solveZoom(const Rocket& rocket, const ZoomOptions& options)
{
    ZoomResult out;
    using Clock = std::chrono::steady_clock;
    auto elapsedSeconds = [](Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double>(end - start).count();
    };

    const std::size_t stages = rocket.stages.size();
    if (stages == 0) {
        return out;
    }

    int fine = options.fine_precision;
    if (fine <= 0) {
        fine = 0;
    }
    if (fine <= 0 || fine > 255) {
        return out;
    }

    int coarse = options.coarse_precision;
    if (coarse <= 0) {
        coarse = std::max(static_cast<int>(stages), fine / 2);
    }
    coarse = std::min(coarse, fine);
    coarse = std::max(coarse, static_cast<int>(stages));

    out.coarse_precision_used = coarse;
    out.fine_precision_used = fine;

    std::size_t top_k = options.top_k;
    if (top_k < 1) {
        top_k = 1;
    }

    FullGridResult coarse_result = solveFullGrid(rocket, coarse, options.verbose, options.threads, top_k, options.full_grid_mode);
    out.coarse_evaluated = coarse_result.combinations;
    out.coarse_best_units = coarse_result.best_units;
    out.coarse_timings = coarse_result.timings;
    if (!std::isfinite(coarse_result.best_mass)) {
        return out;
    }

    const int coarse_radius_steps = std::max(0, options.window_coarse_steps);
    int initial_radius_fine = 0;
    if (coarse_radius_steps > 0) {
        initial_radius_fine = static_cast<int>(std::ceil((static_cast<double>(coarse_radius_steps) * static_cast<double>(fine)) / static_cast<double>(coarse)));
    }

    std::uint64_t budget = options.max_evaluations;
    std::size_t k_effective = top_k;
    std::uint64_t budget_base = budget;
    std::uint64_t budget_extra = 0;

    if (budget != std::numeric_limits<std::uint64_t>::max()) {
        if (budget == 0) {
            budget = 1;
        }
        std::uint64_t k_cap_u64 = budget;
        std::size_t k_cap = (k_cap_u64 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
                                ? std::numeric_limits<std::size_t>::max()
                                : static_cast<std::size_t>(k_cap_u64);
        k_effective = std::min(k_effective, k_cap);
        if (k_effective < 1) {
            k_effective = 1;
        }
        budget_base = budget / static_cast<std::uint64_t>(k_effective);
        budget_extra = budget % static_cast<std::uint64_t>(k_effective);
        if (budget_base == 0) {
            budget_base = 1;
            budget_extra = 0;
        }
    }

    out.best_units.assign(stages, 0);
    out.best_mass = std::numeric_limits<double>::infinity();

    int best_radius_used = 0;
    std::uint64_t refined_total = 0;
    auto refine_start = Clock::now();

    for (std::size_t idx = 0; idx < coarse_result.best_k.size() && idx < k_effective; idx++) {
        const auto& candidate = coarse_result.best_k[idx];
        std::vector<int> center;
        if (!scalePositiveUnitsLargestRemainder(candidate.units.data(), stages, coarse, fine, center)) {
            continue;
        }

        std::uint64_t per_candidate_budget = (budget == std::numeric_limits<std::uint64_t>::max())
                                                 ? std::numeric_limits<std::uint64_t>::max()
                                                 : (budget_base + (idx < budget_extra ? 1ULL : 0ULL));
        RefineResult refined = refineAroundCenter(rocket, fine, center, initial_radius_fine, per_candidate_budget);
        refined_total += refined.evaluated;
        if (!std::isfinite(refined.best_mass)) {
            continue;
        }

        if (betterCandidate(refined.best_mass, refined.best_units, out.best_mass, out.best_units)) {
            out.best_mass = refined.best_mass;
            out.best_units = refined.best_units;
            best_radius_used = refined.used_radius;
        }
    }

    auto refine_end = Clock::now();
    out.refine_seconds = elapsedSeconds(refine_start, refine_end);
    out.radius_used = best_radius_used;
    out.refined_evaluated = refined_total;

    return out;
}

} // namespace optimizer
