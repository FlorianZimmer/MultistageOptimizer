#pragma once

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace util {

inline std::string formatThousands(std::uint64_t value, char separator = ',')
{
    std::string out = std::to_string(value);
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(out.size()) - 3; i > 0; i -= 3) {
        out.insert(static_cast<std::size_t>(i), 1, separator);
    }
    return out;
}

inline std::string formatDoubleThousands(double value, int decimals, char separator = ',')
{
    if (!std::isfinite(value)) {
        return "nan";
    }

    if (decimals < 0) {
        decimals = 0;
    }

    std::uint64_t scale = 1;
    for (int i = 0; i < decimals; i++) {
        scale *= 10ULL;
    }

    bool negative = value < 0.0;
    double abs_value = std::fabs(value);
    double scaled_value = abs_value * static_cast<double>(scale);
    if (scaled_value > static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
        std::ostringstream fallback;
        fallback << std::fixed << std::setprecision(decimals) << value;
        return fallback.str();
    }

    std::uint64_t scaled = static_cast<std::uint64_t>(std::llround(scaled_value));
    std::uint64_t integer = (scale == 0) ? scaled : (scaled / scale);
    std::uint64_t fraction = (scale == 0) ? 0 : (scaled % scale);

    std::ostringstream out;
    if (negative) {
        out << "-";
    }
    out << formatThousands(integer, separator);
    if (decimals > 0) {
        out << "." << std::setw(decimals) << std::setfill('0') << fraction;
    }
    return out.str();
}

} // namespace util

