#pragma once

static json readJson(string path) {
    std::ifstream ifs(path);
    json j = json::parse(ifs);
    return j;
}

static vector<Engine> importEngines(string path) {
    vector<Engine> engineList;
    json file = readJson(path);
    for (int i = 0; i < file["engines"].size(); i++) {
        engineList.push_back(Engine(file["engines"][i]["name"], file["engines"][i]["isp_sl"], file["engines"][i]["isp_vac"], file["engines"][i]["mass"]));
    }
    return engineList;
}

static vector<Engine> find_engines_on_Stage(vector<Engine> engineList, json file) {
    vector<Engine> engines_on_Stage;
    for (int i = 0; i < file["stages"].size(); i++)
    {
        for (int j = 0; j < engineList.size(); j++) {
            if (file["stages"][i]["engine"] == engineList[j].name) {
                engines_on_Stage.push_back(engineList[j]);
            }
        }
    }
    return engines_on_Stage;
}

static vector<Stage> extractStages(vector<Engine> engineList, json file) {
    vector<Stage> all_stages;
    vector<Engine> engines = find_engines_on_Stage(engineList, file);
    for (int i = 0; i < file["stages"].size(); i++) {
        Stage tempStage(0, engines[i], file["stages"][i]["engine_count"], file["stages"][i]["total_mass"], file["stages"][i]["dry_mass"], file["stages"][i]["wet_mass"]);
        all_stages.push_back(tempStage);
    }
    return all_stages;
}

static Rocket importRocket(string path, vector<Engine> engineList) {
    json file = readJson(path);
    double g_body = file["g_body"];
    Rocket rocket(file["g_body"], file["mission_payload"], file["deltaVmission"], extractStages(engineList, file));
    return rocket;
}