#include <string>
#include <vector>

#include "test_framework.hpp"

#include "optimizer/reporting.hpp"

TEST_CASE("buildLowStageDvWarnings: emits warnings for small dv stages")
{
    const std::vector<unsigned char> units = {1, 9};
    auto warnings = optimizer::buildLowStageDvWarnings(units.data(), units.size(), 10, 1000, 200.0);
    REQUIRE(warnings.size() == 1);
    REQUIRE(testfw::stringContains(warnings[0], "Stage 2"));
    REQUIRE(testfw::stringContains(warnings[0], "100.0"));
}

TEST_CASE("buildLowStageDvWarnings: disabled when threshold is 0")
{
    const std::vector<unsigned char> units = {1, 9};
    auto warnings = optimizer::buildLowStageDvWarnings(units.data(), units.size(), 10, 1000, 0.0);
    REQUIRE(warnings.empty());
}

