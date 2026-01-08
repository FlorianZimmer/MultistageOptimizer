#include <string>

#include "test_framework.hpp"

#include "util/budget_config.hpp"

using json = nlohmann::json;

TEST_CASE("util::selectBudgetValue prefers specialized keys")
{
    json config;
    config["maxCombinations"] = 123;
    config["maxCombinationsCpu"] = 456;
    config["maxCombinationsCuda"] = 789;

    const auto* cpu = util::selectBudgetValue(config, util::BudgetBackend::Cpu, "maxCombinations");
    REQUIRE(cpu != nullptr);
    CHECK(cpu->get<int>() == 456);

    const auto* cuda = util::selectBudgetValue(config, util::BudgetBackend::Cuda, "maxCombinations");
    REQUIRE(cuda != nullptr);
    CHECK(cuda->get<int>() == 789);
}

TEST_CASE("util::selectBudgetValue falls back to base key")
{
    json config;
    config["targetSeconds"] = 1.25;

    const auto* cpu = util::selectBudgetValue(config, util::BudgetBackend::Cpu, "targetSeconds");
    REQUIRE(cpu != nullptr);
    CHECK_NEAR(cpu->get<double>(), 1.25, 0.0);
}

TEST_CASE("util::hasAnyBudgetKey detects any budget key")
{
    json config;
    CHECK(!util::hasAnyBudgetKey(config));

    config["maxCombinationsCpu"] = 1;
    CHECK(util::hasAnyBudgetKey(config));
}

TEST_CASE("util::eraseBudgetKeys removes all budget keys")
{
    json config;
    config["maxCombinations"] = 1;
    config["targetSeconds"] = 1.0;
    config["maxCombinationsCpu"] = 1;
    config["targetSecondsCpu"] = 1.0;
    config["maxCombinationsCuda"] = 1;
    config["targetSecondsCuda"] = 1.0;
    config["unrelated"] = 42;

    util::eraseBudgetKeys(config);

    CHECK(!config.contains("maxCombinations"));
    CHECK(!config.contains("targetSeconds"));
    CHECK(!config.contains("maxCombinationsCpu"));
    CHECK(!config.contains("targetSecondsCpu"));
    CHECK(!config.contains("maxCombinationsCuda"));
    CHECK(!config.contains("targetSecondsCuda"));
    CHECK(!util::hasAnyBudgetKey(config));

    CHECK(config.contains("unrelated"));
    CHECK(config["unrelated"].get<int>() == 42);
}

