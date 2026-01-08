#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

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

int main()
{
    const int precision = 5;
    const std::size_t stages = 3;
    const std::uint64_t expected = nCr(precision - 1, static_cast<int>(stages) - 1);
    assert(expected == 6);

    std::vector<unsigned char> distributions(expected * stages, 1);
    std::uint64_t generated = createTuple(distributions.data(), stages, expected, precision);
    assert(generated == expected);
    assert(validateDistributions(distributions.data(), stages, expected, precision, 0));

    std::unordered_set<std::string> seen;
    for (std::uint64_t i = 0; i < expected; i++) {
        auto result = seen.insert(distroKey(distributions.data(), stages, expected, i));
        assert(result.second);
    }

    return 0;
}

