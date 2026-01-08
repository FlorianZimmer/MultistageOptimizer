#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "read_json/json.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "rocket_definition/Stage.hpp"

using json = nlohmann::json;

inline json readJson(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        throw std::runtime_error("Unable to open JSON file: " + path);
    }

    try {
        return json::parse(ifs);
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse JSON file '" + path + "': " + std::string(e.what()));
    }
}

inline std::vector<Engine> importEngines(const std::string& path)
{
    std::vector<Engine> engineList;
    json file = readJson(path);

    const json& engines = file.at("engines");
    engineList.reserve(engines.size());
    for (const auto& engine : engines) {
        engineList.push_back(Engine(engine.at("name").get<std::string>(),
                                    engine.at("isp_sl").get<double>(),
                                    engine.at("isp_vac").get<double>(),
                                    engine.at("mass").get<double>()));
    }
    return engineList;
}

inline std::vector<Stage> extractStages(const std::vector<Engine>& engineList, const json& file)
{
    std::unordered_map<std::string, Engine> enginesByName;
    enginesByName.reserve(engineList.size());
    for (const auto& engine : engineList) {
        enginesByName.emplace(engine.name, engine);
    }

    const json& stages = file.at("stages");
    std::vector<Stage> all_stages;
    all_stages.reserve(stages.size());

    for (const auto& stage : stages) {
        std::string engineName = stage.at("engine").get<std::string>();
        auto it = enginesByName.find(engineName);
        if (it == enginesByName.end()) {
            throw std::runtime_error("Engine '" + engineName + "' not found in engines.json");
        }
        all_stages.push_back(Stage(0.0,
                                   it->second,
                                   stage.at("engine_count").get<int>(),
                                   stage.at("total_mass").get<double>(),
                                   stage.at("dry_mass").get<double>(),
                                   stage.at("wet_mass").get<double>()));
    }

    return all_stages;
}

inline Rocket importRocket(const std::string& path, const std::vector<Engine>& engineList)
{
    json file = readJson(path);
    return Rocket(file.at("g_body").get<double>(),
                  file.at("mission_payload").get<double>(),
                  file.at("deltaVmission").get<int>(),
                  extractStages(engineList, file));
}
