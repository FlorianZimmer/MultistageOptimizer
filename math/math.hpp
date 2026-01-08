#pragma once

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

// Returns value of Binomial Coefficient C(n, k)
// source: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
inline unsigned long long nCr(int n, int k)
{
    if (n < 0 || k < 0 || k > n) {
        return 0;
    }

    // Since C(n, k) = C(n, n-k)
    if (k > n - k) {
        k = n - k;
    }

    std::uint64_t res = 1;
    for (int i = 1; i <= k; i++) {
        std::uint64_t numerator = static_cast<std::uint64_t>(n - (k - i));
        std::uint64_t denominator = static_cast<std::uint64_t>(i);

        std::uint64_t div = std::gcd(numerator, denominator);
        numerator /= div;
        denominator /= div;

        std::uint64_t div2 = std::gcd(res, denominator);
        res /= div2;
        denominator /= div2;

        if (denominator != 1) {
            // Should not happen for binomial coefficients; fall back to clamping.
            return ULLONG_MAX;
        }
        if (numerator != 0 && res > (ULLONG_MAX / numerator)) {
            return ULLONG_MAX;
        }
        res *= numerator;
    }

    return static_cast<unsigned long long>(res);
}

// Fill array with all possible positive-integer tuples ("stars and bars").
// Layout: arr[column * sizeY + row]. Returns the number of generated rows.
inline std::uint64_t createTuple(unsigned char* arr, std::size_t sizeX, std::uint64_t sizeY, int sum)
{
    assert(arr != nullptr);
    assert(sizeX > 0);
    assert(sum >= 0);

    if (sizeX == 0) {
        return 0;
    }

    if (sum < static_cast<int>(sizeX)) {
        return 0;
    }

    if (sum > static_cast<int>(std::numeric_limits<unsigned char>::max())) {
        return 0;
    }

    if (sizeX == 1) {
        if (sizeY < 1) {
            return 0;
        }
        arr[0] = static_cast<unsigned char>(sum);
        return 1;
    }

    const int slots = sum - 1;
    const std::size_t bars_needed = sizeX - 1;
    if (slots < static_cast<int>(bars_needed)) {
        return 0;
    }

    std::vector<int> bars(bars_needed);
    for (std::size_t i = 0; i < bars_needed; i++) {
        bars[i] = static_cast<int>(i) + 1;
    }

    std::uint64_t row = 0;
    for (;;) {
        int previous = 0;
        for (std::size_t col = 0; col < bars_needed; col++) {
            const int b = bars[col];
            const int part = b - previous;
            assert(part >= 1);
            assert(part <= static_cast<int>(std::numeric_limits<unsigned char>::max()));
            arr[col * sizeY + row] = static_cast<unsigned char>(part);
            previous = b;
        }
        const int last = sum - previous;
        assert(last >= 1);
        assert(last <= static_cast<int>(std::numeric_limits<unsigned char>::max()));
        arr[bars_needed * sizeY + row] = static_cast<unsigned char>(last);

#ifndef NDEBUG
        int check_sum = 0;
        for (std::size_t col = 0; col < sizeX; col++) {
            check_sum += arr[col * sizeY + row];
            assert(arr[col * sizeY + row] >= 1);
        }
        assert(check_sum == sum);
#endif

        row++;
        if (row >= sizeY) {
            break;
        }

        int idx = static_cast<int>(bars_needed) - 1;
        while (idx >= 0) {
            const int max_value = slots - (static_cast<int>(bars_needed) - 1 - idx);
            if (bars[static_cast<std::size_t>(idx)] != max_value) {
                break;
            }
            idx--;
        }
        if (idx < 0) {
            break;
        }
        bars[static_cast<std::size_t>(idx)] += 1;
        for (std::size_t j = static_cast<std::size_t>(idx) + 1; j < bars_needed; j++) {
            bars[j] = bars[j - 1] + 1;
        }
    }

    return row;
}

//Addition function that detects overflow
inline unsigned long safeAdd_ULONG(unsigned long lhs, unsigned long rhs)
{
    if (ULONG_MAX - lhs < rhs) {
        return ULONG_MAX;
    }
    return lhs + rhs;
}

inline unsigned long long safeAdd_ULLONG(unsigned long long lhs, unsigned long long rhs)
{
    if (ULLONG_MAX - lhs < rhs) {
        return ULLONG_MAX;
    }
    return lhs + rhs;
}

inline unsigned long long safeMULT_ULLONG(unsigned long long a, unsigned long long b)
{
    if (a != 0 && b > (ULLONG_MAX / a)) {
        return ULLONG_MAX;
    }
    return a * b;
}

inline unsigned long long requiredBytesForPrecision(int precision, int nStages, std::size_t massTypeSize)
{
    (void)massTypeSize;
    unsigned long long combos = nCr(precision - 1, nStages - 1);
    unsigned long long distributions_bytes = safeMULT_ULLONG(combos, static_cast<unsigned long long>(nStages)); // uint8_t
    return distributions_bytes;
}

inline unsigned long long requiredBytesForPrecisionStreaming(int precision, int nStages, std::size_t massTypeSize)
{
    (void)massTypeSize;
    if (precision <= 0 || nStages <= 0) {
        return 0;
    }

    constexpr unsigned long long dbl = static_cast<unsigned long long>(sizeof(double));
    unsigned long long stages = static_cast<unsigned long long>(nStages);

    unsigned long long ratio_entries = safeMULT_ULLONG(stages, static_cast<unsigned long long>(precision + 1));
    unsigned long long ratio_bytes = safeMULT_ULLONG(ratio_entries, dbl);

    unsigned long long per_stage_bytes = safeMULT_ULLONG(stages, 4ULL * dbl); // inert_mass_rest/fract/one_minus/inert_coeff

    return safeAdd_ULLONG(ratio_bytes, per_stage_bytes);
}

//calculate the maximum precision without causing an overflow when creating the array distributions
inline int calcMaxPrecUp(int precision, int nStages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    int i = precision;
    while (i > 1 && requiredBytesForPrecision(i, nStages, massTypeSize) >= maxRAM) {
        i--;
    }
    return i;
}

//calculate the maximum precision without causing an overflow when creating the array distributions
inline int calcMaxPrecDown(int precision, int nStages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    int i = precision;
    while (i > 1 && requiredBytesForPrecision(i, nStages, massTypeSize) >= maxRAM) {
        i--;
    }
    return i;
}

inline int calcMaxPrecUpStreaming(int precision, int nStages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    int i = precision;
    while (i > 1 && requiredBytesForPrecisionStreaming(i, nStages, massTypeSize) >= maxRAM) {
        i--;
    }
    return i;
}

inline int calcMaxPrecDownStreaming(int precision, int nStages, unsigned long long maxRAM, std::size_t massTypeSize)
{
    int i = precision;
    while (i > 1 && requiredBytesForPrecisionStreaming(i, nStages, massTypeSize) >= maxRAM) {
        i--;
    }
    return i;
}

inline double* ArrToPerc(double* newArr, unsigned char* arr, int sizeX, int precision)
{
    for (int i = 0; i < sizeX; i++) {
        newArr[i] = (static_cast<double>(arr[i]) / precision) * 100.0;
    }
    return newArr;
}
