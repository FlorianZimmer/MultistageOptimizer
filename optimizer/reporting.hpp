#pragma once

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace optimizer {

inline std::vector<std::string> buildLowStageDvWarnings(const unsigned char* units,
                                                        std::size_t stages,
                                                        int precision,
                                                        int deltaVmission,
                                                        double minStageDvMps)
{
    std::vector<std::string> warnings;
    if (units == nullptr || stages == 0) {
        return warnings;
    }
    if (precision <= 0) {
        return warnings;
    }
    if (deltaVmission <= 0) {
        return warnings;
    }
    if (!std::isfinite(minStageDvMps) || minStageDvMps <= 0.0) {
        return warnings;
    }

    for (std::size_t i = 0; i < stages; i++) {
        double dv = (static_cast<double>(units[i]) / static_cast<double>(precision)) * static_cast<double>(deltaVmission);
        if (dv >= minStageDvMps) {
            continue;
        }
        double frac = dv / static_cast<double>(deltaVmission);
        int stage_label = static_cast<int>(stages - i);

        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss.precision(1);
        oss << "Stage " << stage_label << " contributes only " << dv << " m/s (" << (frac * 100.0)
            << "%); consider removing/merging this stage.";
        warnings.push_back(oss.str());
    }

    return warnings;
}

} // namespace optimizer

