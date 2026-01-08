#include <cassert>
#include <cmath>

#include "math/massCalc.hpp"

int main()
{
    const double payload = 1000.0;
    const double inert_rest = 250.0;
    const double inert_fract = 0.1;
    const int deltaV = 9400;
    const double g0 = 9.80665;
    const double isp = 350.0;

    double baseline = calcMass(payload, inert_fract, 0.0, deltaV, g0, isp, inert_rest);
    assert(std::isfinite(baseline));
    assert(std::abs(baseline - (payload + inert_rest)) < 1e-9);

    double invalid_inert = calcMass(payload, 1.0, 0.5, deltaV, g0, isp, inert_rest);
    assert(std::isinf(invalid_inert));

    double unreachable = calcMass(payload, 0.9, 1.0, deltaV, g0, isp, inert_rest);
    assert(std::isinf(unreachable));

    double invalid_gravity = calcMass(payload, inert_fract, 0.5, deltaV, 0.0, isp, inert_rest);
    assert(std::isinf(invalid_gravity));

    return 0;
}

