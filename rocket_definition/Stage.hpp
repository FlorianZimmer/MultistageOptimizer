#pragma once

#include "rocket_definition/Engine.hpp"

using namespace std;

class Stage {
public:
    int num_engines;
    Engine engine;
    double tot_mass;
    double wet_mass;
    double dry_mass;
    double prop_mass;  //calc
    double inert_mass_tankonly; //calc
    double inert_mass_rest; //calc
    double inert_mass_fract; //calc
    double payload;
    Stage(double payload, Engine engine, int num_engines, double tot_mass, double dry_mass, double wet_mass) {
        this->payload = payload;

        //can be deleted as well as not necessary info. Weight can be calculated as seen in inert_mass_rest
        this->num_engines = num_engines;
        this->engine = engine;
        this->tot_mass = tot_mass;
        this->wet_mass = wet_mass;
        this->dry_mass = dry_mass;
        prop_mass = wet_mass - dry_mass;
        //calculations are more precise, by knowing what is a scalable inert mass and what not, when increasing or decreasing the fuel in a tank, 
        //that is why the to different variables for inert masses exist, one for calculation of inert mass fraction and the other non scalable, which can include the mass of engine(s), rcs, control unit, etc.
        inert_mass_rest = tot_mass - wet_mass;

        //can be deleted -> dry mass
        inert_mass_tankonly = dry_mass;     //tot_mass - prop_mass - inert_mass_rest;   //inert mass of the fuel tank only
        inert_mass_fract = (double)inert_mass_tankonly / wet_mass;
    };
};

