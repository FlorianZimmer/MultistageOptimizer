#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace optimizer {

inline bool scalePositiveUnitsLargestRemainder(const unsigned char* coarseUnits,
                                               std::size_t stages,
                                               int coarsePrecision,
                                               int finePrecision,
                                               std::vector<int>& outFineUnits)
{
    outFineUnits.clear();
    if (coarseUnits == nullptr || stages == 0) {
        return false;
    }
    if (coarsePrecision <= 0 || finePrecision <= 0) {
        return false;
    }
    if (finePrecision < static_cast<int>(stages)) {
        return false;
    }

    struct Entry {
        int index = 0;
        double fractional = 0.0;
    };

    outFineUnits.assign(stages, 1);
    std::vector<Entry> entries;
    entries.reserve(stages);

    int sum = 0;
    for (std::size_t i = 0; i < stages; i++) {
        const int value = static_cast<int>(coarseUnits[i]);
        if (value <= 0) {
            return false;
        }
        const double scaled = (static_cast<double>(value) * static_cast<double>(finePrecision)) / static_cast<double>(coarsePrecision);
        int base = static_cast<int>(scaled);
        if (base < 1) {
            base = 1;
        }
        outFineUnits[i] = base;
        sum += base;
        entries.push_back(Entry{static_cast<int>(i), scaled - static_cast<double>(static_cast<int>(scaled))});
    }

    auto byFractionDesc = [](const Entry& left, const Entry& right) {
        if (left.fractional != right.fractional) {
            return left.fractional > right.fractional;
        }
        return left.index < right.index;
    };
    auto byFractionAsc = [](const Entry& left, const Entry& right) {
        if (left.fractional != right.fractional) {
            return left.fractional < right.fractional;
        }
        return left.index < right.index;
    };

    int remainder = finePrecision - sum;
    if (remainder > 0) {
        std::sort(entries.begin(), entries.end(), byFractionDesc);
        std::size_t cursor = 0;
        while (remainder > 0) {
            outFineUnits[static_cast<std::size_t>(entries[cursor].index)] += 1;
            remainder--;
            cursor++;
            if (cursor == entries.size()) {
                cursor = 0;
            }
        }
    }
    else if (remainder < 0) {
        std::sort(entries.begin(), entries.end(), byFractionAsc);
        std::size_t cursor = 0;
        while (remainder < 0) {
            int& slot = outFineUnits[static_cast<std::size_t>(entries[cursor].index)];
            if (slot > 1) {
                slot -= 1;
                remainder++;
            }
            cursor++;
            if (cursor == entries.size()) {
                cursor = 0;
            }

            bool any_reducible = false;
            for (int v : outFineUnits) {
                if (v > 1) {
                    any_reducible = true;
                    break;
                }
            }
            if (!any_reducible) {
                return false;
            }
        }
    }

    int check = std::accumulate(outFineUnits.begin(), outFineUnits.end(), 0);
    if (check != finePrecision) {
        return false;
    }
    for (int v : outFineUnits) {
        if (v < 1) {
            return false;
        }
    }
    return true;
}

inline void computeBoundsAroundCenter(const std::vector<int>& centerUnits,
                                      int precision,
                                      int radius,
                                      std::vector<int>& outMinUnits,
                                      std::vector<int>& outMaxUnits)
{
    assert(radius >= 0);
    const std::size_t stages = centerUnits.size();
    outMinUnits.assign(stages, 1);
    outMaxUnits.assign(stages, 1);

    const int global_max = precision - static_cast<int>(stages) + 1;
    for (std::size_t i = 0; i < stages; i++) {
        int center = centerUnits[i];
        outMinUnits[i] = std::max(1, center - radius);
        outMaxUnits[i] = std::min(global_max, center + radius);
    }
}

inline std::uint64_t countBoundedCompositions(const std::vector<int>& minUnits,
                                              const std::vector<int>& maxUnits,
                                              int precision,
                                              std::uint64_t cap = std::numeric_limits<std::uint64_t>::max())
{
    const std::size_t stages = minUnits.size();
    if (stages == 0 || maxUnits.size() != stages) {
        return 0;
    }
    if (precision < static_cast<int>(stages)) {
        return 0;
    }
    for (std::size_t i = 0; i < stages; i++) {
        if (minUnits[i] < 1 || maxUnits[i] < minUnits[i]) {
            return 0;
        }
    }

    std::vector<int> minSuffix(stages + 1, 0);
    std::vector<int> maxSuffix(stages + 1, 0);
    for (std::size_t i = stages; i-- > 0;) {
        minSuffix[i] = minSuffix[i + 1] + minUnits[i];
        maxSuffix[i] = maxSuffix[i + 1] + maxUnits[i];
    }
    if (minSuffix[0] > precision || maxSuffix[0] < precision) {
        return 0;
    }

    std::uint64_t count = 0;
    std::vector<int> stack(stages, 0);
    std::vector<int> value(stages, 0);

    std::size_t idx = 0;
    int remaining = precision;
    value[0] = minUnits[0];

    for (;;) {
        int v = value[idx];
        int remaining_after = remaining - v;
        int min_needed = minSuffix[idx + 1];
        int max_allowed = maxSuffix[idx + 1];

        if (remaining_after >= min_needed && remaining_after <= max_allowed) {
            if (idx + 1 == stages) {
                if (remaining_after == 0) {
                    count++;
                    if (count >= cap) {
                        return cap;
                    }
                }
            }
            else {
                remaining = remaining_after;
                idx++;
                stack[idx] = v;
                value[idx] = minUnits[idx];
                continue;
            }
        }

        for (;;) {
            if (value[idx] < maxUnits[idx]) {
                value[idx]++;
                break;
            }
            if (idx == 0) {
                return count;
            }
            idx--;
            remaining += stack[idx + 1];
        }
    }
}

template <typename Fn>
inline std::uint64_t enumerateBoundedCompositions(const std::vector<int>& minUnits,
                                                  const std::vector<int>& maxUnits,
                                                  int precision,
                                                  Fn&& fn,
                                                  std::uint64_t cap = std::numeric_limits<std::uint64_t>::max())
{
    const std::size_t stages = minUnits.size();
    if (stages == 0 || maxUnits.size() != stages) {
        return 0;
    }

    std::vector<int> minSuffix(stages + 1, 0);
    std::vector<int> maxSuffix(stages + 1, 0);
    for (std::size_t i = stages; i-- > 0;) {
        minSuffix[i] = minSuffix[i + 1] + minUnits[i];
        maxSuffix[i] = maxSuffix[i + 1] + maxUnits[i];
    }
    if (minSuffix[0] > precision || maxSuffix[0] < precision) {
        return 0;
    }

    std::uint64_t generated = 0;
    std::vector<int> units(stages, 0);

    std::vector<int> stack(stages, 0);
    std::vector<int> value(stages, 0);

    std::size_t idx = 0;
    int remaining = precision;
    value[0] = minUnits[0];

    for (;;) {
        int v = value[idx];
        units[idx] = v;
        int remaining_after = remaining - v;
        int min_needed = minSuffix[idx + 1];
        int max_allowed = maxSuffix[idx + 1];

        if (remaining_after >= min_needed && remaining_after <= max_allowed) {
            if (idx + 1 == stages) {
                if (remaining_after == 0) {
                    fn(units);
                    generated++;
                    if (generated >= cap) {
                        return generated;
                    }
                }
            }
            else {
                remaining = remaining_after;
                idx++;
                stack[idx] = v;
                value[idx] = minUnits[idx];
                continue;
            }
        }

        for (;;) {
            if (value[idx] < maxUnits[idx]) {
                value[idx]++;
                break;
            }
            if (idx == 0) {
                return generated;
            }
            idx--;
            remaining += stack[idx + 1];
        }
    }
}

} // namespace optimizer

