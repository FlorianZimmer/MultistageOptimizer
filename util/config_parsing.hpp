#pragma once

#include <cctype>
#include <climits>
#include <cstdlib>
#include <string>

#include "read_json/json.hpp"

inline std::string stripSpacesAndUnderscores(const std::string& value)
{
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '_') {
            continue;
        }
        out.push_back(c);
    }
    return out;
}

inline bool parseByteSizeString(const std::string& raw, unsigned long long& out)
{
    std::string value = stripSpacesAndUnderscores(raw);
    if (value.empty()) {
        return false;
    }

    std::size_t pos = 0;
    while (pos < value.size() && (std::isdigit(static_cast<unsigned char>(value[pos])) || value[pos] == '.')) {
        pos++;
    }

    std::string number_part = value.substr(0, pos);
    std::string suffix = value.substr(pos);
    if (number_part.empty()) {
        return false;
    }

    char* end = nullptr;
    double number = std::strtod(number_part.c_str(), &end);
    if (end == nullptr || *end != '\0' || number < 0.0) {
        return false;
    }

    for (char& c : suffix) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    unsigned long long multiplier = 1;
    if (suffix.empty() || suffix == "b") {
        multiplier = 1;
    }
    else if (suffix == "k" || suffix == "kb" || suffix == "kib") {
        multiplier = 1024ULL;
    }
    else if (suffix == "m" || suffix == "mb" || suffix == "mib") {
        multiplier = 1024ULL * 1024ULL;
    }
    else if (suffix == "g" || suffix == "gb" || suffix == "gib") {
        multiplier = 1024ULL * 1024ULL * 1024ULL;
    }
    else if (suffix == "t" || suffix == "tb" || suffix == "tib") {
        multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL;
    }
    else {
        return false;
    }

    long double bytes = static_cast<long double>(number) * static_cast<long double>(multiplier);
    if (bytes > static_cast<long double>(ULLONG_MAX)) {
        return false;
    }

    out = static_cast<unsigned long long>(bytes);
    return true;
}

inline bool parseMaxRam(const nlohmann::json& value, unsigned long long& out)
{
    if (value.is_number_unsigned()) {
        out = value.get<unsigned long long>();
        return true;
    }
    if (value.is_number_integer()) {
        long long signed_value = value.get<long long>();
        if (signed_value < 0) {
            return false;
        }
        out = static_cast<unsigned long long>(signed_value);
        return true;
    }
    if (value.is_number_float()) {
        double number = value.get<double>();
        if (number < 0.0 || number > static_cast<double>(ULLONG_MAX)) {
            return false;
        }
        out = static_cast<unsigned long long>(number);
        return true;
    }
    if (value.is_string()) {
        return parseByteSizeString(value.get<std::string>(), out);
    }
    return false;
}

