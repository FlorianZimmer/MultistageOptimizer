#pragma once

#include <string>

using namespace std;

class Engine {
public:
    string name;
    double isp_sl;
    double isp_vac;
    double mass;        //kg
    Engine(string name, double isp_sl, double isp_vac, double mass) {
        this->name = name;
        this->isp_sl = isp_sl;
        this->isp_vac = isp_vac;
        this->mass = mass;
    }
    Engine() {
        this->name = "null";
        this->isp_sl = 0;
        this->isp_vac = 0;
        this->mass = 0;
    }
};

