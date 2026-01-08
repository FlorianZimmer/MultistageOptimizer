#include <cstdint>
#include <vector>

#include "test_framework.hpp"

#include "optimizer/bounded_compositions.hpp"

TEST_CASE("scalePositiveUnitsLargestRemainder: preserves sum and positivity")
{
    const std::vector<unsigned char> coarse = {2, 3, 5};
    std::vector<int> fine;
    REQUIRE(optimizer::scalePositiveUnitsLargestRemainder(coarse.data(), coarse.size(), 10, 20, fine));
    REQUIRE(fine.size() == coarse.size());
    REQUIRE(fine[0] == 4);
    REQUIRE(fine[1] == 6);
    REQUIRE(fine[2] == 10);
}

TEST_CASE("scalePositiveUnitsLargestRemainder: distributes remainder deterministically")
{
    const std::vector<unsigned char> coarse = {1, 1, 8};
    std::vector<int> fine;
    REQUIRE(optimizer::scalePositiveUnitsLargestRemainder(coarse.data(), coarse.size(), 10, 15, fine));
    REQUIRE(fine.size() == coarse.size());
    REQUIRE(fine[0] + fine[1] + fine[2] == 15);
    REQUIRE(fine[0] >= 1);
    REQUIRE(fine[1] >= 1);
    REQUIRE(fine[2] >= 1);
    // Both stage 0 and 1 have the same fractional part; tie-break adds to the lower index.
    REQUIRE(fine[0] == 2);
    REQUIRE(fine[1] == 1);
    REQUIRE(fine[2] == 12);
}

TEST_CASE("countBoundedCompositions: matches brute force for small case")
{
    // 3 stages, sum 10
    // bounds:
    //   a in [2,5]
    //   b in [1,6]
    //   c in [1,3]
    std::vector<int> mins = {2, 1, 1};
    std::vector<int> maxs = {5, 6, 3};

    std::uint64_t brute = 0;
    for (int a = mins[0]; a <= maxs[0]; a++) {
        for (int b = mins[1]; b <= maxs[1]; b++) {
            for (int c = mins[2]; c <= maxs[2]; c++) {
                if (a + b + c == 10) {
                    brute++;
                }
            }
        }
    }

    std::uint64_t counted = optimizer::countBoundedCompositions(mins, maxs, 10);
    REQUIRE(counted == brute);
}

TEST_CASE("enumerateBoundedCompositions: generates expected count and respects constraints")
{
    std::vector<int> mins = {1, 2, 1};
    std::vector<int> maxs = {4, 4, 6};
    const int precision = 9;

    std::uint64_t counted = optimizer::countBoundedCompositions(mins, maxs, precision);
    std::uint64_t generated = 0;

    optimizer::enumerateBoundedCompositions(mins, maxs, precision, [&](const std::vector<int>& units) {
        REQUIRE(units.size() == mins.size());
        int sum = 0;
        for (std::size_t i = 0; i < units.size(); i++) {
            REQUIRE(units[i] >= mins[i]);
            REQUIRE(units[i] <= maxs[i]);
            sum += units[i];
        }
        REQUIRE(sum == precision);
        generated++;
    });

    REQUIRE(generated == counted);
}

