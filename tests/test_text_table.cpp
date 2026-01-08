#include "test_framework.hpp"

#include <vector>

#include "util/text_table.hpp"

TEST_CASE("formatTextTable aligns columns and includes separator")
{
    std::vector<std::string> headers = {"metric", "A", "B"};
    std::vector<std::vector<std::string>> rows = {
        {"status", "pass", "fail"},
        {"n", "1", "22"},
    };

    std::string out = util::formatTextTable(headers, rows, 2);

    const std::string expected =
        "metric  A     B   \n"
        "------  ----  ----\n"
        "status  pass  fail\n"
        "n       1     22  \n";

    REQUIRE(out == expected);
}
