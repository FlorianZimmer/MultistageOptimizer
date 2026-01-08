#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

#include "math/math.hpp"

namespace optimizer {

inline std::uint64_t combinationsForPrecision(int precision, int stages)
{
    if (precision <= 0 || stages <= 0) {
        return 0;
    }
    return static_cast<std::uint64_t>(nCr(precision - 1, stages - 1));
}

inline bool fitsFullGridRam(int precision, int stages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    if (precision <= 0 || stages <= 0) {
        return false;
    }
    unsigned long long required = requiredBytesForPrecision(precision, stages, massTypeSize);
    return required < maxRAM;
}

inline bool fitsFullGridRamStreaming(int precision, int stages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    if (precision <= 0 || stages <= 0) {
        return false;
    }
    unsigned long long required = requiredBytesForPrecisionStreaming(precision, stages, massTypeSize);
    return required < maxRAM;
}

inline int choosePrecisionForMaxCombinations(int stages,
                                             std::uint64_t maxCombinations,
                                             unsigned long long maxRAM,
                                             std::size_t massTypeSize,
                                             int maxPrecision = 255)
{
    if (stages <= 0) {
        return 0;
    }
    if (maxCombinations == 0) {
        return 0;
    }
    if (maxPrecision < stages) {
        return 0;
    }

    int low = stages;
    int high = maxPrecision;
    int best = 0;

    while (low <= high) {
        int mid = low + (high - low) / 2;
        std::uint64_t combos = combinationsForPrecision(mid, stages);
        bool combos_ok = (combos != 0 && combos <= maxCombinations);
        bool ram_ok = fitsFullGridRam(mid, stages, maxRAM, massTypeSize);

        if (combos_ok && ram_ok) {
            best = mid;
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }

    return best;
}

inline int choosePrecisionForMaxCombinationsStreaming(int stages,
                                                      std::uint64_t maxCombinations,
                                                      unsigned long long maxRAM,
                                                      std::size_t massTypeSize,
                                                      int maxPrecision = 255)
{
    if (stages <= 0) {
        return 0;
    }
    if (maxCombinations == 0) {
        return 0;
    }
    if (maxPrecision < stages) {
        return 0;
    }

    int low = stages;
    int high = maxPrecision;
    int best = 0;

    while (low <= high) {
        int mid = low + (high - low) / 2;
        std::uint64_t combos = combinationsForPrecision(mid, stages);
        bool combos_ok = (combos != 0 && combos <= maxCombinations);
        bool ram_ok = fitsFullGridRamStreaming(mid, stages, maxRAM, massTypeSize);

        if (combos_ok && ram_ok) {
            best = mid;
            low = mid + 1;
        }
        else {
            high = mid - 1;
        }
    }

    return best;
}

} // namespace optimizer
