#include <climits>

#include "test_framework.hpp"

#include "util/config_parsing.hpp"

TEST_CASE("stripSpacesAndUnderscores: removes whitespace and underscores")
{
    REQUIRE(stripSpacesAndUnderscores(" 16_GB \r\n") == "16GB");
}

TEST_CASE("parseByteSizeString: parses common suffixes (base 1024)")
{
    unsigned long long bytes = 0;

    REQUIRE(parseByteSizeString("1KB", bytes));
    REQUIRE(bytes == 1024ULL);

    REQUIRE(parseByteSizeString("2MiB", bytes));
    REQUIRE(bytes == 2ULL * 1024ULL * 1024ULL);

    REQUIRE(parseByteSizeString("1.5GB", bytes));
    REQUIRE(bytes == 1610612736ULL);
}

TEST_CASE("parseByteSizeString: rejects invalid inputs")
{
    unsigned long long bytes = 0;
    REQUIRE(!parseByteSizeString("", bytes));
    REQUIRE(!parseByteSizeString("   ", bytes));
    REQUIRE(!parseByteSizeString("GB", bytes));
    REQUIRE(!parseByteSizeString("-1GB", bytes));
    REQUIRE(!parseByteSizeString("1XB", bytes));
}

TEST_CASE("parseMaxRam: accepts numbers and strings")
{
    unsigned long long out = 0;

    REQUIRE(parseMaxRam(nlohmann::json(1024), out));
    REQUIRE(out == 1024ULL);

    REQUIRE(parseMaxRam(nlohmann::json(1.5), out));
    REQUIRE(out == 1ULL);

    REQUIRE(!parseMaxRam(nlohmann::json(-1), out));

    REQUIRE(parseMaxRam(nlohmann::json("16GB"), out));
    REQUIRE(out == 16ULL * 1024ULL * 1024ULL * 1024ULL);

    REQUIRE(!parseMaxRam(nlohmann::json(true), out));
}

