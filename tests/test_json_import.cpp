#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "test_framework.hpp"

#include "read_json/read_json.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"

static std::filesystem::path makeTempPath(const char* filename)
{
    namespace fs = std::filesystem;
    fs::path base;
    try {
        base = fs::temp_directory_path();
    }
    catch (...) {
        base = fs::current_path();
    }
    return base / filename;
}

TEST_CASE("readJson: throws for missing file")
{
    const auto missing = makeTempPath("multistageoptimizer_missing.json");
    CHECK_THROWS_WITH(readJson(missing.string()), "Unable to open JSON file");
}

TEST_CASE("readJson: throws for invalid JSON")
{
    namespace fs = std::filesystem;
    const fs::path path = makeTempPath("multistageoptimizer_invalid.json");

    {
        std::ofstream out(path);
        REQUIRE(out.is_open());
        out << R"({"engines":[)";
    }

    CHECK_THROWS_WITH(readJson(path.string()), "Failed to parse JSON file");
    (void)fs::remove(path);
}

TEST_CASE("importEngines: parses engines")
{
    namespace fs = std::filesystem;
    const fs::path path = makeTempPath("multistageoptimizer_engines.json");

    {
        std::ofstream out(path);
        REQUIRE(out.is_open());
        out << R"({"engines":[)"
            << R"({"name":"E1","isp_sl":300,"isp_vac":320,"mass":1},)"
            << R"({"name":"E2","isp_sl":310.5,"isp_vac":330.25,"mass":2.5})"
            << R"(]})";
    }

    const auto engines = importEngines(path.string());
    REQUIRE(engines.size() == 2);
    REQUIRE(engines[0].name == "E1");
    CHECK_NEAR(engines[0].isp_sl, 300.0, 1e-12);
    CHECK_NEAR(engines[0].isp_vac, 320.0, 1e-12);
    CHECK_NEAR(engines[0].mass, 1.0, 1e-12);
    REQUIRE(engines[1].name == "E2");
    CHECK_NEAR(engines[1].isp_sl, 310.5, 1e-12);
    CHECK_NEAR(engines[1].isp_vac, 330.25, 1e-12);
    CHECK_NEAR(engines[1].mass, 2.5, 1e-12);

    (void)fs::remove(path);
}

TEST_CASE("importRocket: throws if stage references unknown engine")
{
    namespace fs = std::filesystem;

    const fs::path engines_path = makeTempPath("multistageoptimizer_tmp_engines.json");
    const fs::path rocket_path = makeTempPath("multistageoptimizer_tmp_rocket_missing_engine.json");

    {
        std::ofstream out(engines_path);
        REQUIRE(out.is_open());
        out << R"({"engines":[{"name":"E1","isp_sl":300,"isp_vac":320,"mass":1}]})";
    }

    {
        std::ofstream out(rocket_path);
        REQUIRE(out.is_open());
        out << R"({"g_body":9.81,"mission_payload":1000,"deltaVmission":1000,"stages":[{"engine":"MISSING","engine_count":1,"total_mass":1,"dry_mass":1,"wet_mass":1}]})";
    }

    const auto engines = importEngines(engines_path.string());
    CHECK_THROWS_WITH(importRocket(rocket_path.string(), engines), "Engine 'MISSING'");

    (void)fs::remove(engines_path);
    (void)fs::remove(rocket_path);
}

TEST_CASE("importRocket: parses rocket and computes stage payload chain")
{
    namespace fs = std::filesystem;

    const fs::path engines_path = makeTempPath("multistageoptimizer_tmp_engines_ok.json");
    const fs::path rocket_path = makeTempPath("multistageoptimizer_tmp_rocket_ok.json");

    {
        std::ofstream out(engines_path);
        REQUIRE(out.is_open());
        out << R"({"engines":[{"name":"E1","isp_sl":300,"isp_vac":320,"mass":1}]})";
    }

    {
        std::ofstream out(rocket_path);
        REQUIRE(out.is_open());
        out << R"({"g_body":9.81,"mission_payload":1000,"deltaVmission":1500,"stages":[)"
            << R"({"engine":"E1","engine_count":1,"total_mass":100,"dry_mass":20,"wet_mass":80},)"
            << R"({"engine":"E1","engine_count":1,"total_mass":300,"dry_mass":60,"wet_mass":240})"
            << R"(]})";
    }

    const auto engines = importEngines(engines_path.string());
    const Rocket rocket = importRocket(rocket_path.string(), engines);

    REQUIRE(rocket.stages.size() == 2);
    CHECK_NEAR(rocket.g_body, 9.81, 1e-12);
    REQUIRE(rocket.deltaV == 1500);

    CHECK_NEAR(rocket.stages[0].payload, 1000.0, 1e-12);
    CHECK_NEAR(rocket.stages[1].payload, rocket.stages[0].tot_mass + rocket.stages[0].payload, 1e-12);

    (void)fs::remove(engines_path);
    (void)fs::remove(rocket_path);
}

