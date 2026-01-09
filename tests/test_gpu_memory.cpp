#include "test_framework.hpp"

#include "util/gpu_memory.hpp"

TEST_CASE("util::sampleProcessGpuMemory returns sane values when available")
{
    util::GpuMemorySample sample = util::sampleProcessGpuMemory();
    if (!sample.available) {
        return;
    }

    CHECK(sample.total_bytes > 0);
    CHECK(sample.used_bytes <= sample.total_bytes);
}

