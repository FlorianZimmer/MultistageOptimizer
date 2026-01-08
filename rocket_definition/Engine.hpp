#pragma once

#include <string>
#include <utility>

class Engine {
public:
    std::string name;
    double isp_sl = 0.0;
    double isp_vac = 0.0;
    double mass = 0.0; // kg

    Engine() = default;
    Engine(std::string name, double isp_sl, double isp_vac, double mass)
        : name(std::move(name)), isp_sl(isp_sl), isp_vac(isp_vac), mass(mass)
    {
    }
};
