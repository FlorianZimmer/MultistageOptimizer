#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "test_framework.hpp"

#include "backend/backend_loader.hpp"
#include "backend/full_grid_problem.hpp"
#include "optimizer/solver.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "rocket_definition/Stage.hpp"

TEST_CASE("backend::parseMode accepts expected values")
{
    backend::Mode mode = backend::Mode::Cpu;
    REQUIRE(backend::parseMode("auto", mode));
    CHECK(mode == backend::Mode::Auto);
    REQUIRE(backend::parseMode("cpu", mode));
    CHECK(mode == backend::Mode::Cpu);
    REQUIRE(backend::parseMode("cuda", mode));
    CHECK(mode == backend::Mode::Cuda);
    REQUIRE(backend::parseMode("rocm", mode));
    CHECK(mode == backend::Mode::Rocm);
}

TEST_CASE("backend::resolve handles auto and unsupported backends")
{
    backend::Request req;
    req.mode = backend::Mode::Auto;
    req.device_id = -1;

    backend::ResolvedBackend resolved;
    std::string error;
    REQUIRE(backend::resolve(req, resolved, error));

    CHECK(resolved.mode == backend::Mode::Cpu || resolved.mode == backend::Mode::Cuda);
    if (resolved.mode == backend::Mode::Cuda) {
        REQUIRE(resolved.plugin.has_value());
        REQUIRE(resolved.plugin->loaded());
    }

    backend::Request rocm_req;
    rocm_req.mode = backend::Mode::Rocm;
    rocm_req.device_id = -1;
    backend::ResolvedBackend rocm_resolved;
    std::string rocm_error;
    CHECK(!backend::resolve(rocm_req, rocm_resolved, rocm_error));
    REQUIRE(!rocm_error.empty());
}

static Rocket makeTestRocket()
{
    Engine e0("E0", 300.0, 450.0, 0.0);
    Engine e1("E1", 300.0, 380.0, 0.0);
    Engine e2("E2", 300.0, 320.0, 0.0);

    std::vector<Stage> stages;
    stages.emplace_back(0.0, e0, 1, 1000.0, 100.0, 900.0);
    stages.emplace_back(0.0, e1, 1, 5000.0, 700.0, 4500.0);
    stages.emplace_back(0.0, e2, 1, 20000.0, 4000.0, 16000.0);

    return Rocket(9.81, 1000.0, 1000, stages);
}

TEST_CASE("CUDA plugin matches CPU streaming full-grid minimum (if available)")
{
    backend::Request req;
    req.mode = backend::Mode::Cuda;
    req.device_id = -1;

    backend::ResolvedBackend resolved;
    std::string error;
    if (!backend::resolve(req, resolved, error) || resolved.mode != backend::Mode::Cuda || !resolved.plugin.has_value()) {
        return; // CUDA plugin not available in this environment.
    }

    Rocket rocket = makeTestRocket();

    const int precision = 10;
    optimizer::FullGridResult cpu = optimizer::solveFullGrid(rocket, precision, false, 1, 1, optimizer::FullGridMode::Streaming);
    REQUIRE(std::isfinite(cpu.best_mass));

    backend::FullGridProblemStorage storage;
    REQUIRE(storage.build(rocket, precision));

    std::vector<std::uint8_t> best_units(storage.problem.stage_count, 0);
    mso_full_grid_solution sol{};
    sol.best_mass = std::numeric_limits<double>::infinity();
    sol.best_index = 0;
    sol.best_units = best_units.data();
    sol.best_units_len = storage.problem.stage_count;

    std::string solve_error;
    mso_status st = resolved.plugin->solveFullGridMin(storage.problem, sol, solve_error);
    REQUIRE_MSG(st == MSO_STATUS_OK, solve_error);
    REQUIRE(std::isfinite(sol.best_mass));

    CHECK_NEAR(sol.best_mass, cpu.best_mass, 1e-12);
    CHECK(sol.best_index == cpu.best_index);

    std::vector<unsigned char> gpu_units(best_units.begin(), best_units.end());
    CHECK(gpu_units == cpu.best_units);

    // Sanity: evaluate the reported distribution on CPU.
    double eval = optimizer::evaluateDistributionMass(rocket, storage.tables, gpu_units.data());
    CHECK_NEAR(eval, sol.best_mass, 1e-12);
}

