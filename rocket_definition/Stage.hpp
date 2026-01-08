#pragma once

#include "rocket_definition/Engine.hpp"
#include <utility>

class Stage {
public:
    int num_engines = 0;
    Engine engine;
    double tot_mass = 0.0;
    double wet_mass = 0.0;
    double dry_mass = 0.0;
    double prop_mass = 0.0;          // calc
    double inert_mass_tankonly = 0.0; // calc
    double inert_mass_rest = 0.0;    // calc
    double inert_mass_fract = 0.0;   // calc
    double payload = 0.0;

    Stage() = default;
    Stage(double payload, Engine engine, int num_engines, double tot_mass, double dry_mass, double wet_mass)
        : num_engines(num_engines),
          engine(std::move(engine)),
          tot_mass(tot_mass),
          wet_mass(wet_mass),
          dry_mass(dry_mass),
          payload(payload)
    {
        prop_mass = wet_mass - dry_mass;
        inert_mass_rest = tot_mass - wet_mass;
        inert_mass_tankonly = dry_mass;
        inert_mass_fract = (wet_mass != 0.0) ? (inert_mass_tankonly / wet_mass) : 0.0;
    }
};
