#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>

inline bool validateDistributions(const unsigned char* distributions,
                                  std::size_t nStages,
                                  std::uint64_t nCombinations,
                                  int precision,
                                  std::uint64_t maxChecks = 100000)
{
    if (distributions == nullptr) {
        return false;
    }
    if (nStages == 0 || nCombinations == 0) {
        return false;
    }
    if (precision <= 0 || precision > 255) {
        return false;
    }

    std::uint64_t checks = nCombinations;
    if (maxChecks > 0) {
        checks = std::min(checks, maxChecks);
    }

    for (std::uint64_t sample = 0; sample < checks; sample++) {
        std::uint64_t i = sample;
        if (checks != nCombinations) {
            if (checks <= 1) {
                i = 0;
            }
            else {
                i = (sample * (nCombinations - 1)) / (checks - 1);
            }
        }

        int sum = 0;
        for (std::size_t stage = 0; stage < nStages; stage++) {
            unsigned char value = distributions[stage * nCombinations + i];
            if (value < 1 || value > precision) {
                return false;
            }
            sum += value;
        }
        if (sum != precision) {
            return false;
        }
    }

    return true;
}
