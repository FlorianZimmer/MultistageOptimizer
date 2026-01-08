#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "test_framework.hpp"

#include "math/math.hpp"
#include "math/validation.hpp"

static std::string distroKey(const unsigned char* distributions,
                            std::size_t stages,
                            std::uint64_t combinations,
                            std::uint64_t index)
{
    std::ostringstream oss;
    for (std::size_t stage = 0; stage < stages; stage++) {
        if (stage != 0) {
            oss << ",";
        }
        oss << static_cast<int>(distributions[stage * combinations + index]);
    }
    return oss.str();
}

TEST_CASE("createTuple: generates expected count and valid distributions")
{
    const int precision = 5;
    const std::size_t stages = 3;
    const std::uint64_t expected = nCr(precision - 1, static_cast<int>(stages) - 1);
    REQUIRE(expected == 6);

    std::vector<unsigned char> distributions(expected * stages, 1);
    std::uint64_t generated = createTuple(distributions.data(), stages, expected, precision);
    REQUIRE(generated == expected);
    REQUIRE(validateDistributions(distributions.data(), stages, expected, precision, 0));

    std::unordered_set<std::string> seen;
    for (std::uint64_t i = 0; i < expected; i++) {
        auto result = seen.insert(distroKey(distributions.data(), stages, expected, i));
        REQUIRE(result.second);
    }
}

TEST_CASE("createTuple: returns 0 when sum too small or too large")
{
    {
        const int sum = 2;
        const std::size_t stages = 3;
        const std::uint64_t combos = 1;
        std::vector<unsigned char> distributions(combos * stages, 0);
        REQUIRE(createTuple(distributions.data(), stages, combos, sum) == 0);
    }

    {
        const int sum = 256;
        const std::size_t stages = 2;
        const std::uint64_t combos = 1;
        std::vector<unsigned char> distributions(combos * stages, 0);
        REQUIRE(createTuple(distributions.data(), stages, combos, sum) == 0);
    }
}

TEST_CASE("createTuple: stops after provided row capacity")
{
    const int precision = 10;
    const std::size_t stages = 3;
    const std::uint64_t expected = nCr(precision - 1, static_cast<int>(stages) - 1);
    REQUIRE(expected == 36);

    const std::uint64_t capacity = 10;
    std::vector<unsigned char> distributions(capacity * stages, 1);
    std::uint64_t generated = createTuple(distributions.data(), stages, capacity, precision);
    REQUIRE(generated == capacity);
    REQUIRE(validateDistributions(distributions.data(), stages, capacity, precision, 0));
}

