#include <vector>

#include "test_framework.hpp"

#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "rocket_definition/Stage.hpp"

TEST_CASE("Stage: computes derived masses")
{
    Engine engine("E1", 300.0, 320.0, 1.0);
    Stage stage(123.0, engine, 2, 1000.0, 100.0, 900.0);

    REQUIRE(stage.num_engines == 2);
    REQUIRE(stage.engine.name == "E1");
    CHECK_NEAR(stage.payload, 123.0, 1e-12);

    CHECK_NEAR(stage.prop_mass, 800.0, 1e-12);
    CHECK_NEAR(stage.inert_mass_rest, 100.0, 1e-12);
    CHECK_NEAR(stage.inert_mass_tankonly, 100.0, 1e-12);
    CHECK_NEAR(stage.inert_mass_fract, 100.0 / 900.0, 1e-12);
}

TEST_CASE("Stage: inert_mass_fract is 0 when wet_mass is 0")
{
    Engine engine("E1", 300.0, 320.0, 1.0);
    Stage stage(0.0, engine, 1, 100.0, 10.0, 0.0);
    CHECK_NEAR(stage.inert_mass_fract, 0.0, 1e-12);
}

TEST_CASE("Rocket: computes payload chain from mission payload and stage total mass")
{
    Engine engine("E1", 300.0, 320.0, 1.0);

    std::vector<Stage> proto;
    proto.emplace_back(0.0, engine, 1, 100.0, 20.0, 80.0);
    proto.emplace_back(0.0, engine, 1, 300.0, 60.0, 240.0);

    Rocket rocket(9.81, 1000.0, 1500, proto);
    REQUIRE(rocket.stages.size() == 2);
    REQUIRE(rocket.deltaV == 1500);
    CHECK_NEAR(rocket.g_body, 9.81, 1e-12);

    CHECK_NEAR(rocket.stages[0].payload, 1000.0, 1e-12);
    CHECK_NEAR(rocket.stages[1].payload, rocket.stages[0].tot_mass + rocket.stages[0].payload, 1e-12);
}

