#include "test_framework.hpp"

#include "util/format_numbers.hpp"

TEST_CASE("formatThousands inserts separators")
{
    REQUIRE(util::formatThousands(0) == "0");
    REQUIRE(util::formatThousands(1) == "1");
    REQUIRE(util::formatThousands(12) == "12");
    REQUIRE(util::formatThousands(123) == "123");
    REQUIRE(util::formatThousands(1234) == "1,234");
    REQUIRE(util::formatThousands(12345) == "12,345");
    REQUIRE(util::formatThousands(123456) == "123,456");
    REQUIRE(util::formatThousands(1234567) == "1,234,567");
}

TEST_CASE("formatDoubleThousands formats decimals and separators")
{
    REQUIRE(util::formatDoubleThousands(6171457.429407, 3) == "6,171,457.429");
    REQUIRE(util::formatDoubleThousands(6171457.429607, 3) == "6,171,457.430");
    REQUIRE(util::formatDoubleThousands(-1234.5, 1) == "-1,234.5");
}

