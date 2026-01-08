#include <cmath>
#include <limits>
#include <vector>

#include "test_framework.hpp"

#include "optimizer/solver.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "rocket_definition/Stage.hpp"

TEST_CASE("solveZoom: wide window matches full-grid optimum mass")
{
    Engine e0("E0", 300.0, 450.0, 0.0);
    Engine e1("E1", 300.0, 380.0, 0.0);
    Engine e2("E2", 300.0, 320.0, 0.0);

    std::vector<Stage> stages;
    stages.emplace_back(0.0, e0, 1, 1000.0, 100.0, 900.0);
    stages.emplace_back(0.0, e1, 1, 5000.0, 700.0, 4500.0);
    stages.emplace_back(0.0, e2, 1, 20000.0, 4000.0, 16000.0);

    Rocket rocket(9.81, 1000.0, 1000, stages);

    const int fine = 10;
    const int coarse = 6;

    optimizer::FullGridResult fine_result = optimizer::solveFullGrid(rocket, fine, false, 1, 1);
    REQUIRE(std::isfinite(fine_result.best_mass));

    optimizer::ZoomOptions opts;
    opts.coarse_precision = coarse;
    opts.fine_precision = fine;
    opts.window_coarse_steps = 100; // large enough to cover full fine grid
    opts.top_k = 2;
    opts.threads = 1;
    opts.verbose = false;
    opts.max_evaluations = std::numeric_limits<std::uint64_t>::max();

    optimizer::ZoomResult zoom_result = optimizer::solveZoom(rocket, opts);
    REQUIRE(std::isfinite(zoom_result.best_mass));
    CHECK_NEAR(zoom_result.best_mass, fine_result.best_mass, 1e-12);
}

TEST_CASE("solveZoom: respects refinement max_evaluations budget")
{
    Engine e0("E0", 300.0, 450.0, 0.0);
    Engine e1("E1", 300.0, 380.0, 0.0);
    Engine e2("E2", 300.0, 320.0, 0.0);

    std::vector<Stage> stages;
    stages.emplace_back(0.0, e0, 1, 1000.0, 100.0, 900.0);
    stages.emplace_back(0.0, e1, 1, 5000.0, 700.0, 4500.0);
    stages.emplace_back(0.0, e2, 1, 20000.0, 4000.0, 16000.0);

    Rocket rocket(9.81, 1000.0, 1000, stages);

    optimizer::ZoomOptions opts;
    opts.coarse_precision = 6;
    opts.fine_precision = 10;
    opts.window_coarse_steps = 100; // would cover full grid without budget
    opts.top_k = 10;
    opts.threads = 1;
    opts.verbose = false;
    opts.max_evaluations = 5;

    optimizer::ZoomResult zoom_result = optimizer::solveZoom(rocket, opts);
    REQUIRE(zoom_result.refined_evaluated <= opts.max_evaluations);
}
