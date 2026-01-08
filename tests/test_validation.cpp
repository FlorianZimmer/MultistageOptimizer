#include <cstdint>
#include <vector>

#include "test_framework.hpp"

#include "math/math.hpp"
#include "math/validation.hpp"

TEST_CASE("validateDistributions: rejects invalid inputs")
{
    REQUIRE(!validateDistributions(nullptr, 1, 1, 10));

    std::vector<unsigned char> dummy(1, 1);
    REQUIRE(!validateDistributions(dummy.data(), 0, 1, 10));
    REQUIRE(!validateDistributions(dummy.data(), 1, 0, 10));
    REQUIRE(!validateDistributions(dummy.data(), 1, 1, 0));
    REQUIRE(!validateDistributions(dummy.data(), 1, 1, 256));
}

TEST_CASE("validateDistributions: accepts valid full table")
{
    const int precision = 5;
    const std::size_t stages = 3;
    const std::uint64_t combos = nCr(precision - 1, static_cast<int>(stages) - 1);
    REQUIRE(combos == 6);

    std::vector<unsigned char> distributions(combos * stages, 0);
    REQUIRE(createTuple(distributions.data(), stages, combos, precision) == combos);
    REQUIRE(validateDistributions(distributions.data(), stages, combos, precision, 0));
}

TEST_CASE("validateDistributions: detects invalid values and invalid sums")
{
    const int precision = 10;
    const std::size_t stages = 2;
    const std::uint64_t combos = nCr(precision - 1, static_cast<int>(stages) - 1);
    REQUIRE(combos == 9);

    std::vector<unsigned char> distributions(combos * stages, 0);
    REQUIRE(createTuple(distributions.data(), stages, combos, precision) == combos);
    REQUIRE(validateDistributions(distributions.data(), stages, combos, precision, 0));

    // Corrupt one entry to 0 (invalid range).
    std::vector<unsigned char> invalid = distributions;
    invalid[0] = 0;
    REQUIRE(!validateDistributions(invalid.data(), stages, combos, precision, 0));

    // Corrupt sum: make (stage0=1, stage1=1) -> sum != precision.
    invalid = distributions;
    invalid[0 * combos + 0] = 1;
    invalid[1 * combos + 0] = 1;
    REQUIRE(!validateDistributions(invalid.data(), stages, combos, precision, 0));
}

TEST_CASE("validateDistributions: sampling can miss a corrupted middle entry")
{
    const int precision = 10;
    const std::size_t stages = 2;
    const std::uint64_t combos = nCr(precision - 1, static_cast<int>(stages) - 1);
    REQUIRE(combos == 9);

    std::vector<unsigned char> distributions(combos * stages, 0);
    REQUIRE(createTuple(distributions.data(), stages, combos, precision) == combos);

    // Corrupt a middle entry only.
    distributions[0 * combos + 4] = 0;
    distributions[1 * combos + 4] = static_cast<unsigned char>(precision);

    // Full validation catches it.
    REQUIRE(!validateDistributions(distributions.data(), stages, combos, precision, 0));

    // Sampling 2 checks only first and last index; the middle corruption is not inspected.
    REQUIRE(validateDistributions(distributions.data(), stages, combos, precision, 2));
}

