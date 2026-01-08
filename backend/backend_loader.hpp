#pragma once

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "backend/mso_backend_api.h"

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

namespace backend {

enum class Mode {
    Auto,
    Cpu,
    Cuda,
    Rocm,
};

struct Request {
    Mode mode = Mode::Auto;
    std::string backend_path; // optional file or directory
    int device_id = -1;
};

inline std::string toLower(std::string_view value)
{
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

inline bool parseMode(std::string_view value, Mode& out)
{
    std::string s = toLower(value);
    if (s == "auto") {
        out = Mode::Auto;
        return true;
    }
    if (s == "cpu") {
        out = Mode::Cpu;
        return true;
    }
    if (s == "cuda" || s == "nvidia") {
        out = Mode::Cuda;
        return true;
    }
    if (s == "rocm" || s == "hip" || s == "amd") {
        out = Mode::Rocm;
        return true;
    }
    return false;
}

inline const char* modeName(Mode mode)
{
    switch (mode) {
    case Mode::Auto:
        return "auto";
    case Mode::Cpu:
        return "cpu";
    case Mode::Cuda:
        return "cuda";
    case Mode::Rocm:
        return "rocm";
    }
    return "unknown";
}

inline std::filesystem::path sharedLibraryFilenameForMode(Mode mode)
{
    const char* base = nullptr;
    switch (mode) {
    case Mode::Cuda:
        base = "mso_backend_cuda";
        break;
    case Mode::Rocm:
        base = "mso_backend_rocm";
        break;
    default:
        base = "";
        break;
    }

    std::filesystem::path name(base);
#if defined(_WIN32)
    name += ".dll";
#elif defined(__APPLE__)
    name += ".dylib";
#else
    name += ".so";
#endif
    return name;
}

inline std::filesystem::path processExecutablePath()
{
#if defined(_WIN32)
    std::string buffer;
    buffer.resize(MAX_PATH);
    for (;;) {
        DWORD n = GetModuleFileNameA(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (n == 0) {
            return {};
        }
        if (n < buffer.size() - 1) {
            buffer.resize(n);
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
    return std::filesystem::path(buffer);
#else
    std::string buffer;
    buffer.resize(1024);
    for (;;) {
        ssize_t n = readlink("/proc/self/exe", buffer.data(), buffer.size());
        if (n < 0) {
            return {};
        }
        if (static_cast<size_t>(n) < buffer.size()) {
            buffer.resize(static_cast<size_t>(n));
            break;
        }
        buffer.resize(buffer.size() * 2);
    }
    return std::filesystem::path(buffer);
#endif
}

inline std::filesystem::path processExecutableDir()
{
    std::filesystem::path exe = processExecutablePath();
    if (exe.empty()) {
        return {};
    }
    return exe.parent_path();
}

inline std::optional<std::string> readEnvVar(std::string_view name)
{
    const char* raw = std::getenv(std::string(name).c_str());
    if (!raw || *raw == '\0') {
        return std::nullopt;
    }
    return std::string(raw);
}

class DynamicLibrary final {
public:
    DynamicLibrary() = default;
    DynamicLibrary(const DynamicLibrary&) = delete;
    DynamicLibrary& operator=(const DynamicLibrary&) = delete;

    DynamicLibrary(DynamicLibrary&& other) noexcept
    {
        *this = std::move(other);
    }

    DynamicLibrary& operator=(DynamicLibrary&& other) noexcept
    {
        if (this == &other) {
            return *this;
        }
        unload();
        handle_ = other.handle_;
        path_ = std::move(other.path_);
        other.handle_ = nullptr;
        return *this;
    }

    ~DynamicLibrary()
    {
        unload();
    }

    bool load(const std::filesystem::path& path, std::string& error)
    {
        unload();
        path_ = path;

#if defined(_WIN32)
        handle_ = LoadLibraryA(path.string().c_str());
        if (!handle_) {
            error = lastErrorString();
            return false;
        }
#else
        handle_ = dlopen(path.string().c_str(), RTLD_NOW);
        if (!handle_) {
            const char* msg = dlerror();
            error = msg ? msg : "dlopen failed";
            return false;
        }
#endif

        return true;
    }

    void unload()
    {
        if (!handle_) {
            return;
        }
#if defined(_WIN32)
        FreeLibrary(handle_);
#else
        dlclose(handle_);
#endif
        handle_ = nullptr;
        path_.clear();
    }

    void* symbol(const char* name, std::string& error) const
    {
        if (!handle_) {
            error = "library not loaded";
            return nullptr;
        }
#if defined(_WIN32)
        FARPROC proc = GetProcAddress(handle_, name);
        if (!proc) {
            error = lastErrorString();
            return nullptr;
        }
        return reinterpret_cast<void*>(proc);
#else
        dlerror();
        void* proc = dlsym(handle_, name);
        const char* msg = dlerror();
        if (msg) {
            error = msg;
            return nullptr;
        }
        return proc;
#endif
    }

    bool loaded() const
    {
        return handle_ != nullptr;
    }

    const std::filesystem::path& path() const
    {
        return path_;
    }

private:
    static std::string lastErrorString()
    {
#if defined(_WIN32)
        DWORD code = GetLastError();
        if (code == 0) {
            return "LoadLibrary/GetProcAddress failed";
        }
        LPSTR buffer = nullptr;
        DWORD size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                                    nullptr,
                                    code,
                                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                                    reinterpret_cast<LPSTR>(&buffer),
                                    0,
                                    nullptr);
        std::string out = (size && buffer) ? std::string(buffer, size) : ("Win32 error " + std::to_string(code));
        if (buffer) {
            LocalFree(buffer);
        }
        while (!out.empty() && (out.back() == '\r' || out.back() == '\n')) {
            out.pop_back();
        }
        return out;
#else
        return "dynamic library error";
#endif
    }

#if defined(_WIN32)
    HMODULE handle_ = nullptr;
#else
    void* handle_ = nullptr;
#endif
    std::filesystem::path path_;
};

class Plugin final {
public:
    Plugin() = default;
    Plugin(const Plugin&) = delete;
    Plugin& operator=(const Plugin&) = delete;

    Plugin(Plugin&& other) noexcept
    {
        *this = std::move(other);
    }

    Plugin& operator=(Plugin&& other) noexcept
    {
        if (this == &other) {
            return *this;
        }
        reset();
        lib_ = std::move(other.lib_);
        api_ = other.api_;
        backend_ = other.backend_;
        kind_ = other.kind_;
        other.api_ = nullptr;
        other.backend_ = nullptr;
        other.kind_ = MSO_BACKEND_KIND_UNKNOWN;
        return *this;
    }

    ~Plugin()
    {
        reset();
    }

    bool loaded() const
    {
        return api_ != nullptr && backend_ != nullptr;
    }

    uint32_t kind() const
    {
        return kind_;
    }

    const char* name() const
    {
        return api_ ? (api_->backend_name ? api_->backend_name : "") : "";
    }

    const std::filesystem::path& path() const
    {
        return lib_.path();
    }

    const char* lastError() const
    {
        if (!api_ || !backend_ || !api_->last_error) {
            return "";
        }
        const char* msg = api_->last_error(backend_);
        return msg ? msg : "";
    }

    mso_status solveFullGridMin(const mso_full_grid_problem& problem, mso_full_grid_solution& solution, std::string& error) const
    {
        if (!loaded() || !api_->solve_full_grid_min) {
            error = "backend not loaded";
            return MSO_STATUS_ERROR;
        }
        mso_status st = api_->solve_full_grid_min(backend_, &problem, &solution);
        if (st != MSO_STATUS_OK) {
            error = lastError();
            if (error.empty()) {
                error = "backend solve_full_grid_min failed";
            }
        }
        return st;
    }

    static std::optional<Plugin> load(const std::filesystem::path& path, int device_id, std::string& error)
    {
        Plugin plugin;

        std::string load_error;
        if (!plugin.lib_.load(path, load_error)) {
            error = "Failed to load backend library '" + path.string() + "': " + load_error;
            return std::nullopt;
        }

        std::string sym_error;
        auto* fn_ptr = reinterpret_cast<mso_backend_get_api_fn>(plugin.lib_.symbol("mso_backend_get_api", sym_error));
        if (!fn_ptr) {
            error = "Backend library '" + path.string() + "' is missing mso_backend_get_api: " + sym_error;
            return std::nullopt;
        }

        const mso_backend_api* api = fn_ptr();
        if (!api) {
            error = "Backend library '" + path.string() + "' returned null API";
            return std::nullopt;
        }
        if (api->api_version != MSO_BACKEND_API_VERSION) {
            error = "Backend library '" + path.string() + "' API version mismatch (expected " +
                    std::to_string(MSO_BACKEND_API_VERSION) + ", got " + std::to_string(api->api_version) + ")";
            return std::nullopt;
        }
        if (!api->create || !api->destroy || !api->solve_full_grid_min) {
            error = "Backend library '" + path.string() + "' has an incomplete API";
            return std::nullopt;
        }

        mso_backend_create_params params{};
        params.device_id = static_cast<int32_t>(device_id);
        params.flags = 0;

        void* backend = nullptr;
        mso_status st = api->create(&params, &backend);
        if (st != MSO_STATUS_OK || !backend) {
            std::string detail;
            if (api->last_error && backend) {
                const char* msg = api->last_error(backend);
                if (msg) {
                    detail = msg;
                }
            }
            if (detail.empty()) {
                detail = "create() failed";
            }
            error = "Backend init failed for '" + path.string() + "': " + detail;
            return std::nullopt;
        }

        plugin.api_ = api;
        plugin.backend_ = backend;
        plugin.kind_ = api->backend_kind;
        return plugin;
    }

    void reset()
    {
        if (api_ && backend_ && api_->destroy) {
            api_->destroy(backend_);
        }
        backend_ = nullptr;
        api_ = nullptr;
        kind_ = MSO_BACKEND_KIND_UNKNOWN;
        lib_.unload();
    }

private:
    DynamicLibrary lib_;
    const mso_backend_api* api_ = nullptr;
    void* backend_ = nullptr;
    uint32_t kind_ = MSO_BACKEND_KIND_UNKNOWN;
};

inline std::vector<std::filesystem::path> defaultSearchRoots(const Request& request)
{
    std::vector<std::filesystem::path> roots;

    if (!request.backend_path.empty()) {
        roots.push_back(std::filesystem::path(request.backend_path));
    }
    if (auto env = readEnvVar("MSO_BACKEND_PATH")) {
        roots.push_back(std::filesystem::path(*env));
    }

    std::filesystem::path exe_dir = processExecutableDir();
    if (!exe_dir.empty()) {
        roots.push_back(exe_dir / "backends");
        roots.push_back(exe_dir);
    }

    return roots;
}

inline bool pathExists(const std::filesystem::path& path)
{
    std::error_code ec;
    return std::filesystem::exists(path, ec);
}

inline std::vector<std::filesystem::path> candidatePluginPaths(const Request& request, Mode mode)
{
    std::vector<std::filesystem::path> out;
    std::filesystem::path filename = sharedLibraryFilenameForMode(mode);
    if (filename.empty()) {
        return out;
    }

    for (const auto& root : defaultSearchRoots(request)) {
        if (root.empty()) {
            continue;
        }
        std::error_code ec;
        bool is_dir = std::filesystem::is_directory(root, ec);
        if (!ec && is_dir) {
            out.push_back(root / filename);
            continue;
        }
        out.push_back(root);
    }

    return out;
}

struct ResolvedBackend {
    Mode mode = Mode::Cpu;
    std::optional<Plugin> plugin;
    std::filesystem::path loaded_path;
    std::string warning;
};

inline bool resolve(const Request& request, ResolvedBackend& out, std::string& error)
{
    out = ResolvedBackend{};

    if (request.mode == Mode::Cpu) {
        out.mode = Mode::Cpu;
        return true;
    }

    auto tryLoadMode = [&](Mode mode, std::string& local_error) -> std::optional<Plugin> {
        for (const auto& candidate : candidatePluginPaths(request, mode)) {
            if (!pathExists(candidate)) {
                continue;
            }
            std::string load_error;
            auto plugin = Plugin::load(candidate, request.device_id, load_error);
            if (plugin) {
                out.loaded_path = candidate;
                return plugin;
            }
            local_error = load_error;
        }
        if (local_error.empty()) {
            local_error = "Backend plugin not found in search paths.";
        }
        return std::nullopt;
    };

    if (request.mode == Mode::Cuda) {
        std::string local_error;
        auto plugin = tryLoadMode(Mode::Cuda, local_error);
        if (!plugin) {
            error = local_error;
            return false;
        }
        out.mode = Mode::Cuda;
        out.plugin = std::move(*plugin);
        return true;
    }

    if (request.mode == Mode::Rocm) {
        error = "ROCm backend not implemented yet.";
        return false;
    }

    // Auto mode: try CUDA then fall back to CPU.
    if (request.mode == Mode::Auto) {
        std::string local_error;
        auto plugin = tryLoadMode(Mode::Cuda, local_error);
        if (plugin) {
            out.mode = Mode::Cuda;
            out.plugin = std::move(*plugin);
            return true;
        }
        out.mode = Mode::Cpu;
        out.warning = local_error;
        return true;
    }

    out.mode = Mode::Cpu;
    return true;
}

} // namespace backend
