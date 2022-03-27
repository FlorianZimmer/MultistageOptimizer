#pragma once

#include <vector>
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Stage.hpp"

using namespace std;

class Rocket {
public:
    vector<Stage> stages;
    int deltaV;
    double g_body;
    Rocket(double g_body, double mission_payload, int deltaV, vector<Stage> stages) //: stages(nStages)
    {
        this->g_body = g_body;
        this->deltaV = deltaV;
        for (int i = 0; i < stages.size(); i++) {
            if (i == 0) {
                Stage newStage = Stage(mission_payload, stages[i].engine, stages[i].num_engines, stages[i].tot_mass, stages[i].dry_mass, stages[i].wet_mass);
                this->stages.push_back(newStage);
            }
            else
            {
                Stage newStage = Stage(this->stages[i - 1].tot_mass + this->stages[i - 1].payload, stages[i].engine, stages[i].num_engines, stages[i].tot_mass, stages[i].dry_mass, stages[i].wet_mass);
                this->stages.push_back(newStage);
            }
        }
    }
};

