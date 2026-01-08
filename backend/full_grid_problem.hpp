#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "backend/mso_backend_api.h"
#include "math/math.hpp"
#include "optimizer/solver.hpp"

namespace backend {

struct FullGridProblemStorage {
    optimizer::RocketEvalTables tables;
    std::vector<double> ratio_tables_flat;
    mso_full_grid_problem problem{};

    bool build(const Rocket& rocket, int precision)
    {
        problem = mso_full_grid_problem{};

        const std::size_t stages = rocket.stages.size();
        if (stages == 0) {
            return false;
        }
        if (precision <= 0 || precision > 255) {
            return false;
        }
        if (stages > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            return false;
        }

        tables = optimizer::buildRocketEvalTables(rocket, precision);
        if (!tables.valid) {
            return false;
        }

        ratio_tables_flat.assign(stages * static_cast<std::size_t>(precision + 1), 0.0);
        for (std::size_t s = 0; s < stages; s++) {
            const auto& src = tables.ratio_tables[s];
            for (int u = 0; u <= precision && static_cast<std::size_t>(u) < src.size(); u++) {
                ratio_tables_flat[s * static_cast<std::size_t>(precision + 1) + static_cast<std::size_t>(u)] = src[static_cast<std::size_t>(u)];
            }
        }

        unsigned long long combos = nCr(precision - 1, static_cast<int>(stages) - 1);

        problem.stage_count = static_cast<std::uint32_t>(stages);
        problem.precision = static_cast<std::uint32_t>(precision);
        problem.combinations = static_cast<std::uint64_t>(combos);
        problem.payload = rocket.stages[0].payload;
        problem.inert_mass_rest = tables.inert_mass_rest.data();
        problem.inert_mass_fract = tables.inert_mass_fract.data();
        problem.one_minus_inert_fract = tables.one_minus_inert_fract.data();
        problem.inert_coeff = tables.inert_coeff.data();
        problem.ratio_tables_flat = ratio_tables_flat.data();
        return true;
    }
};

} // namespace backend

