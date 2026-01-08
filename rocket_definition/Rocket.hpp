#pragma once

#include <vector>
#include "rocket_definition/Stage.hpp"

class Rocket {
public:
    std::vector<Stage> stages;
    int deltaV = 0;
    double g_body = 0.0;

    Rocket() = default;
    Rocket(double g_body, double mission_payload, int deltaV, const std::vector<Stage>& stages) //: stages(nStages)
        : deltaV(deltaV), g_body(g_body)
    {
        this->stages.reserve(stages.size());
        for (std::size_t i = 0; i < stages.size(); i++) {
            double payload = mission_payload;
            if (i > 0) {
                payload = this->stages[i - 1].tot_mass + this->stages[i - 1].payload;
            }
            this->stages.emplace_back(payload,
                                      stages[i].engine,
                                      stages[i].num_engines,
                                      stages[i].tot_mass,
                                      stages[i].dry_mass,
                                      stages[i].wet_mass);
        }
    }
};
