#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "test_framework.hpp"

#include "math/massCalc.hpp"
#include "math/math.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "rocket_definition/Stage.hpp"

TEST_CASE("massRatio: 0 dV fraction yields ratio 1")
{
    CHECK_NEAR(massRatio(0.0, 9400, 9.80665, 350.0), 1.0, 1e-12);
}

TEST_CASE("massRatio: increases with deltaV fraction")
{
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;

    const double r1 = massRatio(0.1, deltaV, g0, isp);
    const double r2 = massRatio(0.2, deltaV, g0, isp);
    REQUIRE(std::isfinite(r1));
    REQUIRE(std::isfinite(r2));
    REQUIRE(r2 > r1);
}

TEST_CASE("calcMassFromRatio: ratio=1 returns payload + inert_rest")
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;

    const double one_minus_inert = 1.0 - inert_fract;
    const double inert_coeff = inert_fract / one_minus_inert;

    const double total = calcMassFromRatio(payload, inert_fract, one_minus_inert, inert_coeff, 1.0, inert_rest);
    REQUIRE(std::isfinite(total));
    CHECK_NEAR(total, payload + inert_rest, 1e-12);
}

TEST_CASE("calcMassFromRatio: returns inf when denominator is non-positive")
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    const double payload = 1000.0;
    const double inert_rest = 0.0;
    const double inert_fract = 0.5;

    const double one_minus_inert = 1.0 - inert_fract;
    const double inert_coeff = inert_fract / one_minus_inert;

    // denom = 1 - inert_fract * ratio = 0
    CHECK(calcMassFromRatio(payload, inert_fract, one_minus_inert, inert_coeff, 2.0, inert_rest) == inf);
}

TEST_CASE("calcMass: baseline (0 dV fraction) returns payload + inert_rest")
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;

    const double total = calcMass(payload, inert_fract, 0.0, deltaV, g0, isp, inert_rest);
    REQUIRE(std::isfinite(total));
    CHECK_NEAR(total, payload + inert_rest, 1e-12);
}

TEST_CASE("calcMass: rejects invalid inputs")
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;
    constexpr double inf = std::numeric_limits<double>::infinity();

    CHECK(calcMass(-1.0, inert_fract, 0.1, deltaV, g0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, inert_fract, 0.1, deltaV, g0, isp, -1.0) == inf);
    CHECK(calcMass(payload, 0.0, 0.1, deltaV, g0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, 1.0, 0.1, deltaV, g0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, inert_fract, -0.1, deltaV, g0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, inert_fract, 0.1, -1, g0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, inert_fract, 0.1, deltaV, 0.0, isp, inert_rest) == inf);
    CHECK(calcMass(payload, inert_fract, 0.1, deltaV, g0, 0.0, inert_rest) == inf);
}

TEST_CASE("calcMass: returns inf when inertMassFract * ratio >= 1")
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;

    const double total = calcMass(payload, 0.9, 1.0, deltaV, g0, isp, inert_rest);
    CHECK(total == inf);
}

TEST_CASE("calcMass: matches calcMassFromRatio")
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;
    const double dv_fraction = 0.4;

    const double ratio = massRatio(dv_fraction, deltaV, g0, isp);
    REQUIRE(std::isfinite(ratio));

    const double one_minus_inert = 1.0 - inert_fract;
    const double inert_coeff = inert_fract / one_minus_inert;

    const double expected = calcMassFromRatio(payload, inert_fract, one_minus_inert, inert_coeff, ratio, inert_rest);
    const double actual = calcMass(payload, inert_fract, dv_fraction, deltaV, g0, isp, inert_rest);

    REQUIRE(std::isfinite(expected));
    REQUIRE(std::isfinite(actual));
    CHECK_NEAR(actual, expected, 1e-12);
}

TEST_CASE("calcMass: increases with deltaV fraction")
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;

    const double m1 = calcMass(payload, inert_fract, 0.2, deltaV, g0, isp, inert_rest);
    const double m2 = calcMass(payload, inert_fract, 0.4, deltaV, g0, isp, inert_rest);
    REQUIRE(std::isfinite(m1));
    REQUIRE(std::isfinite(m2));
    REQUIRE(m2 > m1);
}

TEST_CASE("distMass: marks all entries invalid for bad precision")
{
    Engine engine("E", 300.0, 350.0, 1.0);
    std::vector<Stage> proto;
    proto.emplace_back(0.0, engine, 1, 1000.0, 100.0, 900.0);
    Rocket rocket(9.80665, 1000.0, 1000, proto);

    constexpr std::uint64_t nCombinations = 1;
    std::vector<unsigned char> distributions(nCombinations * rocket.stages.size(), 1);
    std::vector<double> masses(nCombinations, 0.0);

    distMass(false, 0, nCombinations, nCombinations, rocket, masses.data(), distributions.data(), 0);
    REQUIRE(std::isinf(masses[0]));
}

TEST_CASE("distMass: 1-stage rocket matches calcMass at full deltaV")
{
    Engine engine("E", 300.0, 350.0, 1.0);
    std::vector<Stage> proto;
    proto.emplace_back(0.0, engine, 1, 1000.0, 100.0, 900.0);

    const double g_body = 9.80665;
    const double mission_payload = 1000.0;
    const int deltaV = 1000;
    Rocket rocket(g_body, mission_payload, deltaV, proto);

    const int precision = 10;
    const std::uint64_t nCombinations = nCr(precision - 1, static_cast<int>(rocket.stages.size()) - 1);
    REQUIRE(nCombinations == 1);

    std::vector<unsigned char> distributions(nCombinations * rocket.stages.size(), 0);
    std::uint64_t generated = createTuple(distributions.data(), rocket.stages.size(), nCombinations, precision);
    REQUIRE(generated == 1);

    std::vector<double> masses(nCombinations, 0.0);
    distMass(false, 0, nCombinations, nCombinations, rocket, masses.data(), distributions.data(), precision);

    const Stage& stage = rocket.stages[0];
    const double expected = calcMass(stage.payload,
                                    stage.inert_mass_fract,
                                    1.0,
                                    rocket.deltaV,
                                    rocket.g_body,
                                    stage.engine.isp_vac,
                                    stage.inert_mass_rest);

    REQUIRE(std::isfinite(masses[0]));
    REQUIRE(std::isfinite(expected));
    CHECK_NEAR(masses[0], expected, 1e-12);
}

