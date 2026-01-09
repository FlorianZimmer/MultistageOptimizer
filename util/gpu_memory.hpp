#pragma once

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#ifdef near
#undef near
#endif
#ifdef far
#undef far
#endif
#else
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace util {

struct GpuMemorySample {
    std::uint64_t used_bytes = 0;
    std::uint64_t total_bytes = 0;
    bool available = false;
};

namespace detail {

struct NvmlApi final {
    using nvmlReturn_t = int;
    struct nvmlDevice_st;
    using nvmlDevice_t = nvmlDevice_st*;

    struct nvmlMemory_t {
        unsigned long long total = 0;
        unsigned long long free = 0;
        unsigned long long used = 0;
    };

    struct nvmlProcessInfo_t {
        unsigned int pid = 0;
        unsigned long long usedGpuMemory = 0;
    };

    static constexpr nvmlReturn_t NVML_SUCCESS = 0;
    static constexpr nvmlReturn_t NVML_ERROR_INSUFFICIENT_SIZE = 7;
    static constexpr unsigned long long NVML_VALUE_NOT_AVAILABLE = 0xFFFFFFFFFFFFFFFFull;

#if defined(_WIN32)
    HMODULE lib = nullptr;
#else
    void* lib = nullptr;
#endif
    bool initialized = false;

    nvmlReturn_t (*nvmlInit)() = nullptr;
    nvmlReturn_t (*nvmlDeviceGetCount)(unsigned int* count) = nullptr;
    nvmlReturn_t (*nvmlDeviceGetHandleByIndex)(unsigned int index, nvmlDevice_t* device) = nullptr;
    nvmlReturn_t (*nvmlDeviceGetMemoryInfo)(nvmlDevice_t device, nvmlMemory_t* memory) = nullptr;
    nvmlReturn_t (*nvmlDeviceGetComputeRunningProcesses)(nvmlDevice_t device, unsigned int* infoCount, nvmlProcessInfo_t* infos) = nullptr;
    nvmlReturn_t (*nvmlDeviceGetGraphicsRunningProcesses)(nvmlDevice_t device, unsigned int* infoCount, nvmlProcessInfo_t* infos) = nullptr;

    void init()
    {
#if defined(_WIN32)
        lib = LoadLibraryA("nvml.dll");
#else
        lib = dlopen("libnvidia-ml.so.1", RTLD_NOW);
        if (!lib) {
            lib = dlopen("libnvidia-ml.so", RTLD_NOW);
        }
#endif
        if (!lib) {
            return;
        }

        auto loadSym = [&](const char* name) -> void* {
#if defined(_WIN32)
            return reinterpret_cast<void*>(GetProcAddress(lib, name));
#else
            return dlsym(lib, name);
#endif
        };

        nvmlInit = reinterpret_cast<nvmlReturn_t (*)()>(loadSym("nvmlInit_v2"));
        if (!nvmlInit) {
            nvmlInit = reinterpret_cast<nvmlReturn_t (*)()>(loadSym("nvmlInit"));
        }

        nvmlDeviceGetCount = reinterpret_cast<nvmlReturn_t (*)(unsigned int*)>(loadSym("nvmlDeviceGetCount_v2"));
        if (!nvmlDeviceGetCount) {
            nvmlDeviceGetCount = reinterpret_cast<nvmlReturn_t (*)(unsigned int*)>(loadSym("nvmlDeviceGetCount"));
        }

        nvmlDeviceGetHandleByIndex =
            reinterpret_cast<nvmlReturn_t (*)(unsigned int, nvmlDevice_t*)>(loadSym("nvmlDeviceGetHandleByIndex_v2"));
        if (!nvmlDeviceGetHandleByIndex) {
            nvmlDeviceGetHandleByIndex =
                reinterpret_cast<nvmlReturn_t (*)(unsigned int, nvmlDevice_t*)>(loadSym("nvmlDeviceGetHandleByIndex"));
        }

        nvmlDeviceGetMemoryInfo = reinterpret_cast<nvmlReturn_t (*)(nvmlDevice_t, nvmlMemory_t*)>(loadSym("nvmlDeviceGetMemoryInfo"));

        nvmlDeviceGetComputeRunningProcesses =
            reinterpret_cast<nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*)>(loadSym("nvmlDeviceGetComputeRunningProcesses_v2"));
        if (!nvmlDeviceGetComputeRunningProcesses) {
            nvmlDeviceGetComputeRunningProcesses =
                reinterpret_cast<nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*)>(loadSym("nvmlDeviceGetComputeRunningProcesses"));
        }

        nvmlDeviceGetGraphicsRunningProcesses =
            reinterpret_cast<nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*)>(loadSym("nvmlDeviceGetGraphicsRunningProcesses_v2"));
        if (!nvmlDeviceGetGraphicsRunningProcesses) {
            nvmlDeviceGetGraphicsRunningProcesses =
                reinterpret_cast<nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*)>(loadSym("nvmlDeviceGetGraphicsRunningProcesses"));
        }

        if (!nvmlInit || !nvmlDeviceGetCount || !nvmlDeviceGetHandleByIndex) {
            return;
        }
        if (nvmlInit() != NVML_SUCCESS) {
            return;
        }

        initialized = true;
    }

    std::uint64_t queryProcessUsedBytes(nvmlDevice_t device,
                                       unsigned int pid,
                                       nvmlReturn_t (*fn)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*)) const
    {
        if (!initialized || !fn) {
            return 0;
        }

        thread_local std::vector<nvmlProcessInfo_t> processes(64);
        unsigned int count = static_cast<unsigned int>(processes.size());
        nvmlReturn_t st = fn(device, &count, processes.data());
        if (st == NVML_ERROR_INSUFFICIENT_SIZE) {
            processes.resize(count);
            count = static_cast<unsigned int>(processes.size());
            st = fn(device, &count, processes.data());
        }
        if (st != NVML_SUCCESS) {
            return 0;
        }

        std::uint64_t used = 0;
        for (unsigned int i = 0; i < count; i++) {
            if (processes[i].pid != pid) {
                continue;
            }
            if (processes[i].usedGpuMemory == NVML_VALUE_NOT_AVAILABLE) {
                continue;
            }
            used += static_cast<std::uint64_t>(processes[i].usedGpuMemory);
        }
        return used;
    }
};

inline NvmlApi& nvml()
{
    static NvmlApi api;
    static std::once_flag once;
    std::call_once(once, [&] { api.init(); });
    return api;
}

} // namespace detail

inline GpuMemorySample sampleProcessGpuMemory()
{
    auto& api = detail::nvml();
    if (!api.initialized || !api.nvmlDeviceGetMemoryInfo) {
        return {};
    }

    unsigned int device_count = 0;
    if (api.nvmlDeviceGetCount(&device_count) != detail::NvmlApi::NVML_SUCCESS) {
        return {};
    }
    if (device_count == 0) {
        return {};
    }

    std::uint64_t total_bytes = 0;
    std::uint64_t used_bytes = 0;

    for (unsigned int device_index = 0; device_index < device_count; device_index++) {
        detail::NvmlApi::nvmlDevice_t device = nullptr;
        if (api.nvmlDeviceGetHandleByIndex(device_index, &device) != detail::NvmlApi::NVML_SUCCESS) {
            continue;
        }

        detail::NvmlApi::nvmlMemory_t memory{};
        if (api.nvmlDeviceGetMemoryInfo(device, &memory) == detail::NvmlApi::NVML_SUCCESS) {
            total_bytes += static_cast<std::uint64_t>(memory.total);
            used_bytes += static_cast<std::uint64_t>(memory.used);
        }
    }

    GpuMemorySample sample;
    sample.available = true;
    sample.used_bytes = used_bytes;
    sample.total_bytes = total_bytes;
    return sample;
}

} // namespace util
