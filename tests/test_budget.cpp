#include <cstdint>
#include <limits>

#include "test_framework.hpp"

#include "optimizer/budget.hpp"

TEST_CASE("choosePrecisionForMaxCombinations: picks largest precision that fits")
{
    const int stages = 3;
    const std::uint64_t maxCombos = 10;
    const unsigned long long maxRAM = std::numeric_limits<unsigned long long>::max();
    const std::size_t massTypeSize = sizeof(double);

    int chosen = optimizer::choosePrecisionForMaxCombinations(stages, maxCombos, maxRAM, massTypeSize, 255);
    REQUIRE(chosen == 6); // nCr(5,2)=10, nCr(6,2)=15 would exceed budget.
    REQUIRE(optimizer::combinationsForPrecision(chosen, stages) <= maxCombos);
    REQUIRE(optimizer::combinationsForPrecision(chosen + 1, stages) > maxCombos);
}

TEST_CASE("choosePrecisionForMaxCombinations: respects maxPrecision cap")
{
    const int stages = 3;
    const std::uint64_t maxCombos = 1000000;
    const unsigned long long maxRAM = std::numeric_limits<unsigned long long>::max();
    const std::size_t massTypeSize = sizeof(double);

    int chosen = optimizer::choosePrecisionForMaxCombinations(stages, maxCombos, maxRAM, massTypeSize, 10);
    REQUIRE(chosen <= 10);
    REQUIRE(chosen >= stages);
}

