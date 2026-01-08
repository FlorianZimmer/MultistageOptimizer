#pragma once

#include <string>
#include <string_view>

#include "read_json/json.hpp"

namespace util {

enum class BudgetBackend {
    Cpu,
    Cuda,
};

inline std::string budgetKey(std::string_view base, BudgetBackend backend)
{
    std::string key(base);
    if (backend == BudgetBackend::Cuda) {
        key += "Cuda";
    }
    else {
        key += "Cpu";
    }
    return key;
}

inline const nlohmann::json* selectBudgetValue(const nlohmann::json& config, BudgetBackend backend, const char* base_key)
{
    std::string specialized = budgetKey(base_key, backend);
    if (config.contains(specialized)) {
        return &config.at(specialized);
    }
    if (config.contains(base_key)) {
        return &config.at(base_key);
    }
    return nullptr;
}

inline bool hasAnyBudgetKey(const nlohmann::json& config)
{
    return config.contains("maxCombinations") || config.contains("targetSeconds") ||
           config.contains("maxCombinationsCpu") || config.contains("targetSecondsCpu") ||
           config.contains("maxCombinationsCuda") || config.contains("targetSecondsCuda");
}

inline void eraseBudgetKeys(nlohmann::json& config)
{
    config.erase("maxCombinations");
    config.erase("targetSeconds");
    config.erase("maxCombinationsCpu");
    config.erase("targetSecondsCpu");
    config.erase("maxCombinationsCuda");
    config.erase("targetSecondsCuda");
}

} // namespace util
