#pragma once

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

#include "rocket_definition/Rocket.hpp"

inline double massRatio(double deltaVfraction, int deltaVmission, double gravity, double isp)
{
    double exponent = (deltaVfraction * static_cast<double>(deltaVmission)) / (gravity * isp);
    return std::exp(exponent);
}

inline double calcMassFromRatio(double payloadMass,
                                double inertMassFract,
                                double oneMinusInertMassFract,
                                double inertMassCoeff,
                                double ratio,
                                double inert_mass_rest)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    double denom = 1.0 - (inertMassFract * ratio);
    if (denom <= 0.0) {
        return inf;
    }

    double effective_payload = payloadMass + inert_mass_rest;
    double propMass = (effective_payload * ((ratio - 1.0) * oneMinusInertMassFract)) / denom;
    if (propMass < 0.0) {
        return inf;
    }

    double inertMass = inertMassCoeff * propMass + inert_mass_rest;
    return payloadMass + propMass + inertMass;
}

inline double calcMass(double payloadMass,
                       double inertMassFract,
                       double deltaVfraction,
                       int deltaVmission,
                       double gravity,
                       double isp,
                       double inert_mass_rest)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    if (!std::isfinite(payloadMass) || payloadMass < 0.0) {
        return inf;
    }
    if (!std::isfinite(inert_mass_rest) || inert_mass_rest < 0.0) {
        return inf;
    }
    if (!std::isfinite(inertMassFract) || inertMassFract <= 0.0 || inertMassFract >= 1.0) {
        return inf;
    }
    if (!std::isfinite(deltaVfraction) || deltaVfraction < 0.0) {
        return inf;
    }
    if (deltaVmission < 0) {
        return inf;
    }
    if (!std::isfinite(gravity) || gravity <= 0.0) {
        return inf;
    }
    if (!std::isfinite(isp) || isp <= 0.0) {
        return inf;
    }

    double ratio = massRatio(deltaVfraction, deltaVmission, gravity, isp);
    if (!std::isfinite(ratio) || ratio <= 0.0) {
        return inf;
    }
    if (inertMassFract * ratio >= 1.0) {
        return inf;
    }

    double oneMinusInert = 1.0 - inertMassFract;
    double inertCoeff = inertMassFract / oneMinusInert;
    double total = calcMassFromRatio(payloadMass, inertMassFract, oneMinusInert, inertCoeff, ratio, inert_mass_rest);
    if (!std::isfinite(total) || total < 0.0) {
        return inf;
    }
    return total;
}

//calculate total mass of rocket for all distributions
static void distMass(bool verbose,
                     std::uint64_t start,
                     std::uint64_t end,
                     std::uint64_t nCombinations,
                     const Rocket& rocket,
                     double* massDistributions,
                     const unsigned char* distributions,
                     int precision)
{
    double min = std::numeric_limits<double>::infinity();

    //getting all constant vars of the rocket here as to not repeatedly call them over an over for every iteration in for-loop -> speed ups program by about 30%
    std::size_t rocket_size = rocket.stages.size();
    int deltaV = rocket.deltaV;
    double g_body = rocket.g_body;
    auto markInvalid = [&]() {
        for (std::uint64_t i = start; i < end; i++) {
            massDistributions[static_cast<std::size_t>(i)] = std::numeric_limits<double>::infinity();
        }
    };
    if (precision <= 0) {
        markInvalid();
        return;
    }
    if (!std::isfinite(g_body) || g_body <= 0.0) {
        markInvalid();
        return;
    }
    if (deltaV < 0) {
        markInvalid();
        return;
    }

    std::vector<double> inert_mass_rest(rocket_size);
    std::vector<double> inert_mass_fract(rocket_size);
    std::vector<double> one_minus_inert_fract(rocket_size);
    std::vector<double> inert_coeff(rocket_size);
    std::vector<std::vector<double>> ratio_tables(rocket_size);
    for (std::size_t j = 0; j < rocket_size; j++) {
        inert_mass_rest[j] = rocket.stages[j].inert_mass_rest;
        inert_mass_fract[j] = rocket.stages[j].inert_mass_fract;
        double isp_vac = rocket.stages[j].engine.isp_vac;
        if (!std::isfinite(inert_mass_rest[j]) || inert_mass_rest[j] < 0.0) {
            markInvalid();
            return;
        }
        if (!std::isfinite(inert_mass_fract[j]) || inert_mass_fract[j] <= 0.0 || inert_mass_fract[j] >= 1.0) {
            markInvalid();
            return;
        }
        if (!std::isfinite(isp_vac) || isp_vac <= 0.0) {
            markInvalid();
            return;
        }

        one_minus_inert_fract[j] = 1.0 - inert_mass_fract[j];
        inert_coeff[j] = inert_mass_fract[j] / one_minus_inert_fract[j];

        ratio_tables[j].resize(static_cast<std::size_t>(precision) + 1);
        double coeff = static_cast<double>(deltaV) / (static_cast<double>(precision) * g_body * isp_vac);
        for (int v = 0; v <= precision; v++) {
            double ratio = std::exp(static_cast<double>(v) * coeff);
            ratio_tables[j][static_cast<std::size_t>(v)] = ratio;
        }
    }

    for (std::uint64_t i = start; i < end; i++) {
        //if (i == 258100748){
        //        __debugbreak();
        //}
        for (std::size_t j = 0; j < rocket_size; j++) {
            std::size_t offset = j * static_cast<std::size_t>(nCombinations) + static_cast<std::size_t>(i);
            unsigned int dv_units = distributions[offset];
            if (dv_units > static_cast<unsigned int>(precision)) {
                massDistributions[static_cast<std::size_t>(i)] = std::numeric_limits<double>::infinity();
                continue;
            }
            double ratio = ratio_tables[j][dv_units];
            if (j == 0) {
                massDistributions[static_cast<std::size_t>(i)] =
                    calcMassFromRatio(rocket.stages[j].payload,
                                      inert_mass_fract[j],
                                      one_minus_inert_fract[j],
                                      inert_coeff[j],
                                      ratio,
                                      inert_mass_rest[j]); // uppermost stage only has mission payload as payload
            }
            //else if (j == rocket.stages.size() - 1) {
            //    massDistributions[i] += calcMass(rocket.stages[j].payload, rocket.stages[j].inert_mass_fract, distributions[j * nCombinations + i] / (double)precision, rocket.deltaV, rocket.g_body, rocket.stages[j].engine.isp_sl, rocket.stages[j]);    //first stage assumes surface level isp
            //}
            else {
                massDistributions[static_cast<std::size_t>(i)] =
                    calcMassFromRatio(massDistributions[static_cast<std::size_t>(i)],
                                      inert_mass_fract[j],
                                      one_minus_inert_fract[j],
                                      inert_coeff[j],
                                      ratio,
                                      inert_mass_rest[j]);
            }
        }
        //output for debbugging
        //only outputs into console, when new best mass is found
        //only coherent with single-threaded execution
        if (verbose == true) {
            if (massDistributions[static_cast<std::size_t>(i)] < min) {
                min = massDistributions[static_cast<std::size_t>(i)];
                std::cout << "Distribution number: " << i << "\t" << massDistributions[static_cast<std::size_t>(i)] << "\n";
                for (std::size_t j = 0; j < rocket_size; j++) {
                    std::size_t offset = j * static_cast<std::size_t>(nCombinations) + static_cast<std::size_t>(i);
                    std::cout << "Stage " << rocket_size - j << ": " << static_cast<int>(distributions[offset]) << "\t";
                }
                std::cout << "\n\n";
            }
        }
    }
    /*std::cout << "thread " << (int)((double)((double)end / (double)nCombinations) * 24) << " closed\n";*/
}
