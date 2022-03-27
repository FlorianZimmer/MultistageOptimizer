#pragma once

#include <iostream>

#include "rocket_definition/Rocket.hpp"
#include "math/math.hpp"

using namespace std;

//mass ratio
static double massRatio(double deltaVfraction, int deltaVmission, double gravity, double isp) {
    return std::exp((deltaVfraction * deltaVmission) / (gravity * isp));
}

static unsigned long calcMass(unsigned long payloadMass, float inertMassFract, float deltaVfraction, int deltaVmission, float gravity, float isp, double inert_mass_rest) {
    double R = massRatio(deltaVfraction, deltaVmission, gravity, isp);
    unsigned long mass = 0;

    //when inertMassFract * R >= 1 than the targeted deltav is phisically not achievable with the given effiency of the engine and the inert mass fraction of the tank
    if (inertMassFract * R >= 1) {
        mass = ULONG_MAX;
    }
    else {
        //safe addition as overflows occure here regularly, if overflow is detected the result is ULONG_MAX
        double propMass = ((safeAdd_ULONG(payloadMass, inert_mass_rest) * ((R - 1) * (1 - inertMassFract)))) / (1 - (inertMassFract * R));   //adding non-scalable inert mass as payload
        unsigned long inertMass = safeAdd_ULONG(((inertMassFract / (1 - inertMassFract)) * propMass), inert_mass_rest);
        mass = safeAdd_ULONG(safeAdd_ULONG(propMass, inertMass), payloadMass);
    }

    return mass;
}

//calculate total mass of rocket for all distributions
static void distMass(bool verbose, int start, int end, long long nCombinations, Rocket rocket, unsigned long* massDistributions, unsigned char* distributions, int precision) {
    int min = INT_MAX;

    //getting all constant vars of the rocket here as to not repeatedly call them over an over for every iteration in for-loop -> speed ups program by about 30%
    int rocket_size = rocket.stages.size();
    int deltaV = rocket.deltaV;
    double g_body = rocket.g_body;
    double* inert_mass_rest = new double[rocket_size];
    double* inert_mass_fract = new double[rocket_size];
    double* isp_vac = new double[rocket_size];
    for (int j = 0; j < rocket_size; j++) {
        inert_mass_rest[j] = rocket.stages[j].inert_mass_rest;
        inert_mass_fract[j] = rocket.stages[j].inert_mass_fract;
        isp_vac[j] = rocket.stages[j].engine.isp_vac;
    }

    for (int i = start; i < end; i++) {
        //if (i == 258100748){
        //        __debugbreak();
        //}
        for (int j = 0; j < rocket_size; j++) {
            if (j == 0) {
                massDistributions[i] = calcMass(rocket.stages[j].payload, inert_mass_fract[j], distributions[j * nCombinations + i] / (double)precision, deltaV, g_body, isp_vac[j], inert_mass_rest[j]);    //uppermost stage only has mission payload as payload
            }
            //else if (j == rocket.stages.size() - 1) {
            //    massDistributions[i] += calcMass(rocket.stages[j].payload, rocket.stages[j].inert_mass_fract, distributions[j * nCombinations + i] / (double)precision, rocket.deltaV, rocket.g_body, rocket.stages[j].engine.isp_sl, rocket.stages[j]);    //first stage assumes surface level isp
            //}
            else {
                massDistributions[i] = calcMass(massDistributions[i], inert_mass_fract[j], distributions[j * nCombinations + i] / (double)precision, deltaV, g_body, isp_vac[j], inert_mass_rest[j]);
            }
        }
        //output for debbugging
        //only outputs into console, when new best mass is found
        //only coherent with single-threaded execution
        if (verbose == true) {
            if (massDistributions[i] < min) {
                min = massDistributions[i];
                std::cout << "Distribution number: " << i << "\t" << massDistributions[i] << "\n";
                for (int j = 0; j < rocket.stages.size(); j++) {
                    std::cout << "Stage " << rocket_size - j << ": " << (int)distributions[j * nCombinations + i] << "\t";
                }
                std::cout << "\n\n";
            }
        }
    }
    /*std::cout << "thread " << (int)((double)((double)end / (double)nCombinations) * 24) << " closed\n";*/
}