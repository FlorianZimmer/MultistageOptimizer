#pragma once

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace util {

inline std::string formatTextTable(const std::vector<std::string>& headers,
                                   const std::vector<std::vector<std::string>>& rows,
                                   std::size_t padding_spaces = 2)
{
    if (headers.empty()) {
        return std::string{};
    }

    std::vector<std::size_t> widths(headers.size(), 0);
    for (std::size_t col = 0; col < headers.size(); col++) {
        widths[col] = headers[col].size();
    }
    for (const auto& row : rows) {
        for (std::size_t col = 0; col < headers.size(); col++) {
            const std::string_view cell = (col < row.size()) ? std::string_view(row[col]) : std::string_view{};
            widths[col] = std::max(widths[col], cell.size());
        }
    }

    std::string padding(padding_spaces, ' ');

    auto writeRow = [&](std::ostringstream& out, const std::vector<std::string>& cells) {
        for (std::size_t col = 0; col < headers.size(); col++) {
            if (col > 0) {
                out << padding;
            }
            const std::string& text = (col < cells.size()) ? cells[col] : std::string{};
            out << std::left << std::setw(static_cast<int>(widths[col])) << text;
        }
        out << "\n";
    };

    std::ostringstream out;
    writeRow(out, headers);

    std::vector<std::string> separator;
    separator.reserve(headers.size());
    for (std::size_t col = 0; col < headers.size(); col++) {
        separator.push_back(std::string(widths[col], '-'));
    }
    writeRow(out, separator);

    for (const auto& row : rows) {
        writeRow(out, row);
    }

    return out.str();
}

} // namespace util
