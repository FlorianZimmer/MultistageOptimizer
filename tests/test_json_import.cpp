#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

#include "read_json/read_json.hpp"

int main()
{
    namespace fs = std::filesystem;

    fs::path engines_path = fs::path("tests") / "tmp_engines.json";
    fs::path rocket_path = fs::path("tests") / "tmp_rocket.json";

    {
        std::ofstream out(engines_path);
        out << R"({"engines":[{"name":"E1","isp_sl":300,"isp_vac":320,"mass":1}]})";
    }

    {
        std::ofstream out(rocket_path);
        out << R"({"g_body":9.81,"mission_payload":1000,"deltaVmission":1000,"stages":[{"engine":"MISSING","engine_count":1,"total_mass":1,"dry_mass":1,"wet_mass":1}]})";
    }

    bool threw = false;
    try {
        auto engines = importEngines(engines_path.string());
        auto rocket = importRocket(rocket_path.string(), engines);
        (void)rocket;
    }
    catch (const std::exception& e) {
        threw = true;
        std::string message = e.what();
        assert(message.find("Engine 'MISSING'") != std::string::npos);
    }

    fs::remove(engines_path);
    fs::remove(rocket_path);

    assert(threw);
    return 0;
}

