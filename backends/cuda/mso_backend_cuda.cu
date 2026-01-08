#include <cuda_runtime.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "backend/mso_backend_api.h"

namespace {

constexpr int kMaxStages = 64;

inline const char* cudaErrorString(cudaError_t err)
{
    const char* msg = cudaGetErrorString(err);
    return msg ? msg : "unknown CUDA error";
}

struct CudaBackend {
    int device_id = -1;
    std::string last_error;

    int cached_stage_count = 0;
    int cached_precision = 0;

    double* d_inert_mass_rest = nullptr;
    double* d_inert_mass_fract = nullptr;
    double* d_one_minus_inert_fract = nullptr;
    double* d_inert_coeff = nullptr;
    double* d_ratio_tables = nullptr;
    std::uint64_t* d_binom = nullptr;
    int binom_stride = 0;

    double* d_block_best_mass = nullptr;
    std::uint64_t* d_block_best_index = nullptr;
    std::size_t block_capacity = 0;
};

inline void setError(CudaBackend* backend, std::string message)
{
    if (backend) {
        backend->last_error = std::move(message);
    }
}

inline void clearError(CudaBackend* backend)
{
    if (backend) {
        backend->last_error.clear();
    }
}

inline std::uint64_t clampAddU64(std::uint64_t a, std::uint64_t b)
{
    constexpr std::uint64_t maxv = std::numeric_limits<std::uint64_t>::max();
    if (maxv - a < b) {
        return maxv;
    }
    return a + b;
}

static std::vector<std::uint64_t> buildBinomTable(int precision, int stage_count, int& stride_out)
{
    // Table for C(n, k) with n in [0, precision-1] and k in [0, stage_count-1].
    const int n_max = std::max(0, precision - 1);
    const int k_max = std::max(0, stage_count - 1);

    const int stride = std::max(1, stage_count);
    stride_out = stride;

    std::vector<std::uint64_t> binom(static_cast<std::size_t>((n_max + 1) * stride), 0);
    for (int n = 0; n <= n_max; n++) {
        binom[static_cast<std::size_t>(n * stride + 0)] = 1;
        const int upto = std::min(n, k_max);
        for (int k = 1; k <= upto; k++) {
            if (k == n) {
                binom[static_cast<std::size_t>(n * stride + k)] = 1;
                continue;
            }
            std::uint64_t left = binom[static_cast<std::size_t>((n - 1) * stride + (k - 1))];
            std::uint64_t right = binom[static_cast<std::size_t>((n - 1) * stride + k)];
            binom[static_cast<std::size_t>(n * stride + k)] = clampAddU64(left, right);
        }
    }
    return binom;
}

inline std::uint64_t lookupBinom(const std::vector<std::uint64_t>& binom, int stride, int n, int k)
{
    if (n < 0 || k < 0) {
        return 0;
    }
    std::size_t idx = static_cast<std::size_t>(n * stride + k);
    if (idx >= binom.size()) {
        return 0;
    }
    return binom[idx];
}

static bool unrankCompositionPositive(std::uint64_t index,
                                     int stage_count,
                                     int precision,
                                     const std::vector<std::uint64_t>& binom,
                                     int stride,
                                     std::uint8_t* out_units)
{
    if (!out_units || stage_count <= 0 || precision <= 0) {
        return false;
    }
    if (stage_count > kMaxStages) {
        return false;
    }
    if (precision < stage_count) {
        return false;
    }

    if (stage_count == 1) {
        out_units[0] = static_cast<std::uint8_t>(precision);
        return true;
    }

    int remaining_sum = precision;
    std::uint64_t rank = index;
    for (int s = 0; s < stage_count; s++) {
        int remaining_stages = stage_count - s;
        int unit = 1;
        if (remaining_stages == 1) {
            unit = remaining_sum;
        }
        else {
            int max_unit = remaining_sum - (remaining_stages - 1);
            for (int u = 1; u <= max_unit; u++) {
                int rest_sum = remaining_sum - u;
                int rest_stages = remaining_stages - 1;
                std::uint64_t count = 0;
                if (rest_stages == 1) {
                    count = 1;
                }
                else {
                    count = lookupBinom(binom, stride, rest_sum - 1, rest_stages - 1);
                }
                if (rank >= count) {
                    rank -= count;
                    continue;
                }
                unit = u;
                break;
            }
        }
        if (unit < 1 || unit > 255) {
            return false;
        }
        out_units[s] = static_cast<std::uint8_t>(unit);
        remaining_sum -= unit;
    }
    return remaining_sum == 0;
}

__device__ __forceinline__ bool finite(double x)
{
    return isfinite(x);
}

__device__ __forceinline__ double infd()
{
    return __longlong_as_double(0x7ff0000000000000ULL);
}

__device__ __forceinline__ std::uint64_t binomDevice(const std::uint64_t* binom, int stride, int n, int k)
{
    if (n < 0 || k < 0) {
        return 0;
    }
    return binom[n * stride + k];
}

__device__ __forceinline__ double calcMassFromRatioDevice(double payloadMass,
                                                          double inertMassFract,
                                                          double oneMinusInertMassFract,
                                                          double inertMassCoeff,
                                                          double ratio,
                                                          double inert_mass_rest)
{
    const double denom = 1.0 - (inertMassFract * ratio);
    if (denom <= 0.0) {
        return infd();
    }

    const double effective_payload = payloadMass + inert_mass_rest;
    const double propMass = (effective_payload * ((ratio - 1.0) * oneMinusInertMassFract)) / denom;
    if (propMass < 0.0) {
        return infd();
    }

    const double inertMass = inertMassCoeff * propMass + inert_mass_rest;
    return payloadMass + propMass + inertMass;
}

struct Candidate {
    double mass;
    std::uint64_t index;
};

__device__ __forceinline__ Candidate better(Candidate a, Candidate b)
{
    if (a.mass < b.mass) {
        return a;
    }
    if (a.mass > b.mass) {
        return b;
    }
    return (a.index < b.index) ? a : b;
}

__global__ void fullGridMinKernel(std::uint64_t n_combinations,
                                  int stage_count,
                                  int precision,
                                  double payload,
                                  const double* inert_mass_rest,
                                  const double* inert_mass_fract,
                                  const double* one_minus_inert_fract,
                                  const double* inert_coeff,
                                  const double* ratio_tables,
                                  const std::uint64_t* binom,
                                  int binom_stride,
                                  double* block_best_mass,
                                  std::uint64_t* block_best_index)
{
    std::uint64_t global = static_cast<std::uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    Candidate cand;
    cand.index = global;
    cand.mass = infd();

    if (global < n_combinations) {
        if (stage_count == 1) {
            int unit = precision;
            double ratio = ratio_tables[unit];
            cand.mass = calcMassFromRatioDevice(payload,
                                                inert_mass_fract[0],
                                                one_minus_inert_fract[0],
                                                inert_coeff[0],
                                                ratio,
                                                inert_mass_rest[0]);
        }
        else {
            std::uint64_t rank = global;
            int remaining_sum = precision;
            double current = payload;

            for (int s = 0; s < stage_count; s++) {
                int remaining_stages = stage_count - s;
                int unit = 1;
                if (remaining_stages == 1) {
                    unit = remaining_sum;
                }
                else {
                    int max_unit = remaining_sum - (remaining_stages - 1);
                    for (int u = 1; u <= max_unit; u++) {
                        int rest_sum = remaining_sum - u;
                        int rest_stages = remaining_stages - 1;
                        std::uint64_t count = 0;
                        if (rest_stages == 1) {
                            count = 1;
                        }
                        else {
                            count = binomDevice(binom, binom_stride, rest_sum - 1, rest_stages - 1);
                        }
                        if (rank >= count) {
                            rank -= count;
                            continue;
                        }
                        unit = u;
                        break;
                    }
                }

                remaining_sum -= unit;

                double ratio = ratio_tables[s * (precision + 1) + unit];
                current = calcMassFromRatioDevice(current,
                                                  inert_mass_fract[s],
                                                  one_minus_inert_fract[s],
                                                  inert_coeff[s],
                                                  ratio,
                                                  inert_mass_rest[s]);
                if (!finite(current)) {
                    current = infd();
                    break;
                }
            }

            cand.mass = current;
        }
    }
    else {
        cand.index = 0xFFFFFFFFFFFFFFFFULL;
    }

    extern __shared__ unsigned char scratch[];
    auto* shared_mass = reinterpret_cast<double*>(scratch);
    auto* shared_index = reinterpret_cast<std::uint64_t*>(shared_mass + blockDim.x);

    shared_mass[threadIdx.x] = cand.mass;
    shared_index[threadIdx.x] = cand.index;
    __syncthreads();

    for (unsigned int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            Candidate a{shared_mass[threadIdx.x], shared_index[threadIdx.x]};
            Candidate b{shared_mass[threadIdx.x + offset], shared_index[threadIdx.x + offset]};
            Candidate best = better(a, b);
            shared_mass[threadIdx.x] = best.mass;
            shared_index[threadIdx.x] = best.index;
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        block_best_mass[blockIdx.x] = shared_mass[0];
        block_best_index[blockIdx.x] = shared_index[0];
    }
}

static void destroyDeviceBuffers(CudaBackend& backend)
{
    cudaFree(backend.d_inert_mass_rest);
    cudaFree(backend.d_inert_mass_fract);
    cudaFree(backend.d_one_minus_inert_fract);
    cudaFree(backend.d_inert_coeff);
    cudaFree(backend.d_ratio_tables);
    cudaFree(backend.d_binom);
    cudaFree(backend.d_block_best_mass);
    cudaFree(backend.d_block_best_index);

    backend.d_inert_mass_rest = nullptr;
    backend.d_inert_mass_fract = nullptr;
    backend.d_one_minus_inert_fract = nullptr;
    backend.d_inert_coeff = nullptr;
    backend.d_ratio_tables = nullptr;
    backend.d_binom = nullptr;
    backend.binom_stride = 0;
    backend.d_block_best_mass = nullptr;
    backend.d_block_best_index = nullptr;
    backend.block_capacity = 0;

    backend.cached_stage_count = 0;
    backend.cached_precision = 0;
}

static bool ensureDeviceProblemBuffers(CudaBackend& backend,
                                       int stage_count,
                                       int precision,
                                       const std::vector<std::uint64_t>& binom,
                                       int binom_stride,
                                       std::size_t blocks,
                                       std::string& error)
{
    if (stage_count <= 0 || precision <= 0) {
        error = "Invalid stage_count/precision";
        return false;
    }

    bool need_realloc = (backend.cached_stage_count != stage_count) || (backend.cached_precision != precision) ||
                        (backend.binom_stride != binom_stride);
    if (need_realloc) {
        destroyDeviceBuffers(backend);
    }

    auto ensureAlloc = [&](void** ptr, std::size_t bytes) -> bool {
        if (*ptr) {
            return true;
        }
        cudaError_t st = cudaMalloc(ptr, bytes);
        if (st != cudaSuccess) {
            error = cudaErrorString(st);
            return false;
        }
        return true;
    };

    std::size_t stage_bytes = static_cast<std::size_t>(stage_count) * sizeof(double);
    std::size_t ratio_bytes = static_cast<std::size_t>(stage_count) * static_cast<std::size_t>(precision + 1) * sizeof(double);
    std::size_t binom_bytes = binom.size() * sizeof(std::uint64_t);

    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_inert_mass_rest), stage_bytes)) {
        return false;
    }
    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_inert_mass_fract), stage_bytes)) {
        return false;
    }
    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_one_minus_inert_fract), stage_bytes)) {
        return false;
    }
    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_inert_coeff), stage_bytes)) {
        return false;
    }
    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_ratio_tables), ratio_bytes)) {
        return false;
    }
    if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_binom), binom_bytes)) {
        return false;
    }

    if (blocks > backend.block_capacity) {
        cudaFree(backend.d_block_best_mass);
        cudaFree(backend.d_block_best_index);
        backend.d_block_best_mass = nullptr;
        backend.d_block_best_index = nullptr;
        backend.block_capacity = 0;

        std::size_t out_mass_bytes = blocks * sizeof(double);
        std::size_t out_index_bytes = blocks * sizeof(std::uint64_t);
        if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_block_best_mass), out_mass_bytes)) {
            return false;
        }
        if (!ensureAlloc(reinterpret_cast<void**>(&backend.d_block_best_index), out_index_bytes)) {
            return false;
        }
        backend.block_capacity = blocks;
    }

    backend.cached_stage_count = stage_count;
    backend.cached_precision = precision;
    backend.binom_stride = binom_stride;
    return true;
}

static bool copyProblemToDevice(CudaBackend& backend,
                                const mso_full_grid_problem& problem,
                                const std::vector<std::uint64_t>& binom,
                                std::string& error)
{
    int stage_count = static_cast<int>(problem.stage_count);
    int precision = static_cast<int>(problem.precision);

    std::size_t stage_bytes = static_cast<std::size_t>(stage_count) * sizeof(double);
    std::size_t ratio_bytes = static_cast<std::size_t>(stage_count) * static_cast<std::size_t>(precision + 1) * sizeof(double);

    cudaError_t st = cudaMemcpy(backend.d_inert_mass_rest, problem.inert_mass_rest, stage_bytes, cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaMemcpy(backend.d_inert_mass_fract, problem.inert_mass_fract, stage_bytes, cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaMemcpy(backend.d_one_minus_inert_fract, problem.one_minus_inert_fract, stage_bytes, cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaMemcpy(backend.d_inert_coeff, problem.inert_coeff, stage_bytes, cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaMemcpy(backend.d_ratio_tables, problem.ratio_tables_flat, ratio_bytes, cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }

    st = cudaMemcpy(backend.d_binom, binom.data(), binom.size() * sizeof(std::uint64_t), cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }

    return true;
}

static bool launchSolve(CudaBackend& backend,
                        const mso_full_grid_problem& problem,
                        const std::vector<std::uint64_t>& binom,
                        int binom_stride,
                        mso_full_grid_solution& solution,
                        std::string& error)
{
    const std::uint64_t n_combinations = problem.combinations;
    const int stage_count = static_cast<int>(problem.stage_count);
    const int precision = static_cast<int>(problem.precision);

    if (n_combinations == 0) {
        error = "Invalid combinations count.";
        return false;
    }

    constexpr int threads = 256;
    std::size_t blocks = static_cast<std::size_t>((n_combinations + threads - 1) / threads);

    if (!ensureDeviceProblemBuffers(backend, stage_count, precision, binom, binom_stride, blocks, error)) {
        return false;
    }
    if (!copyProblemToDevice(backend, problem, binom, error)) {
        return false;
    }

    std::size_t shared_bytes = static_cast<std::size_t>(threads) * (sizeof(double) + sizeof(std::uint64_t));

    fullGridMinKernel<<<static_cast<unsigned int>(blocks), threads, shared_bytes>>>(
        n_combinations,
        stage_count,
        precision,
        problem.payload,
        backend.d_inert_mass_rest,
        backend.d_inert_mass_fract,
        backend.d_one_minus_inert_fract,
        backend.d_inert_coeff,
        backend.d_ratio_tables,
        backend.d_binom,
        binom_stride,
        backend.d_block_best_mass,
        backend.d_block_best_index);

    cudaError_t st = cudaGetLastError();
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaDeviceSynchronize();
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }

    std::vector<double> h_mass(blocks, std::numeric_limits<double>::infinity());
    std::vector<std::uint64_t> h_index(blocks, std::numeric_limits<std::uint64_t>::max());
    st = cudaMemcpy(h_mass.data(), backend.d_block_best_mass, blocks * sizeof(double), cudaMemcpyDeviceToHost);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }
    st = cudaMemcpy(h_index.data(), backend.d_block_best_index, blocks * sizeof(std::uint64_t), cudaMemcpyDeviceToHost);
    if (st != cudaSuccess) {
        error = cudaErrorString(st);
        return false;
    }

    double best_mass = std::numeric_limits<double>::infinity();
    std::uint64_t best_index = std::numeric_limits<std::uint64_t>::max();
    for (std::size_t i = 0; i < blocks; i++) {
        double m = h_mass[i];
        std::uint64_t idx = h_index[i];
        if (m < best_mass) {
            best_mass = m;
            best_index = idx;
            continue;
        }
        if (m == best_mass && idx < best_index) {
            best_index = idx;
        }
    }

    solution.best_mass = best_mass;
    solution.best_index = best_index;

    if (solution.best_units && solution.best_units_len >= problem.stage_count) {
        if (!unrankCompositionPositive(best_index, stage_count, precision, binom, binom_stride, solution.best_units)) {
            error = "Failed to unrank best_index to units.";
            return false;
        }
    }

    return true;
}

static std::uint64_t computeCombinations(int precision, int stage_count, const std::vector<std::uint64_t>& binom, int stride)
{
    if (stage_count <= 0 || precision <= 0) {
        return 0;
    }
    if (stage_count == 1) {
        return 1;
    }
    // Positive compositions: C(precision-1, stage_count-1)
    return lookupBinom(binom, stride, precision - 1, stage_count - 1);
}

static mso_status createBackend(const mso_backend_create_params* params, void** out_backend)
{
    if (!out_backend) {
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    *out_backend = nullptr;

    int device_id = -1;
    if (params) {
        if (params->flags != 0) {
            return MSO_STATUS_INVALID_ARGUMENT;
        }
        device_id = params->device_id;
    }

    int device_count = 0;
    cudaError_t st = cudaGetDeviceCount(&device_count);
    if (st != cudaSuccess) {
        return MSO_STATUS_NO_DEVICE;
    }
    if (device_count <= 0) {
        return MSO_STATUS_NO_DEVICE;
    }
    if (device_id < 0) {
        device_id = 0;
    }
    if (device_id >= device_count) {
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    st = cudaSetDevice(device_id);
    if (st != cudaSuccess) {
        return MSO_STATUS_ERROR;
    }
    st = cudaFree(nullptr);
    if (st != cudaSuccess) {
        return MSO_STATUS_ERROR;
    }

    auto* backend = new (std::nothrow) CudaBackend();
    if (!backend) {
        return MSO_STATUS_OUT_OF_MEMORY;
    }
    backend->device_id = device_id;
    *out_backend = backend;
    return MSO_STATUS_OK;
}

static void destroyBackend(void* backend_ptr)
{
    auto* backend = static_cast<CudaBackend*>(backend_ptr);
    if (!backend) {
        return;
    }
    destroyDeviceBuffers(*backend);
    delete backend;
}

static const char* lastError(void* backend_ptr)
{
    auto* backend = static_cast<CudaBackend*>(backend_ptr);
    if (!backend) {
        return "";
    }
    return backend->last_error.c_str();
}

static mso_status solveFullGridMin(void* backend_ptr,
                                  const mso_full_grid_problem* problem,
                                  mso_full_grid_solution* solution)
{
    auto* backend = static_cast<CudaBackend*>(backend_ptr);
    if (!backend || !problem || !solution) {
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    clearError(backend);

    if (problem->stage_count == 0 || problem->precision == 0) {
        setError(backend, "stage_count and precision must be positive");
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    if (problem->stage_count > static_cast<std::uint32_t>(kMaxStages)) {
        setError(backend, "stage_count too large");
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    if (problem->precision > 255) {
        setError(backend, "precision must be <= 255");
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    if (problem->precision < problem->stage_count) {
        setError(backend, "precision must be >= stage_count");
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    if (!problem->inert_mass_rest || !problem->inert_mass_fract || !problem->one_minus_inert_fract || !problem->inert_coeff ||
        !problem->ratio_tables_flat) {
        setError(backend, "problem tables must be non-null");
        return MSO_STATUS_INVALID_ARGUMENT;
    }
    if (solution->best_units && solution->best_units_len < problem->stage_count) {
        setError(backend, "best_units_len must be >= stage_count");
        return MSO_STATUS_INVALID_ARGUMENT;
    }

    cudaError_t st = cudaSetDevice(backend->device_id);
    if (st != cudaSuccess) {
        setError(backend, cudaErrorString(st));
        return MSO_STATUS_ERROR;
    }

    int binom_stride = 0;
    std::vector<std::uint64_t> binom = buildBinomTable(static_cast<int>(problem->precision),
                                                       static_cast<int>(problem->stage_count),
                                                       binom_stride);
    std::uint64_t combinations = problem->combinations;
    if (combinations == 0) {
        combinations = computeCombinations(static_cast<int>(problem->precision), static_cast<int>(problem->stage_count), binom, binom_stride);
    }
    if (combinations == 0 || combinations == std::numeric_limits<std::uint64_t>::max()) {
        setError(backend, "Unsupported combinations count");
        return MSO_STATUS_UNSUPPORTED;
    }

    mso_full_grid_problem patched = *problem;
    patched.combinations = combinations;

    std::string error;
    if (!launchSolve(*backend, patched, binom, binom_stride, *solution, error)) {
        setError(backend, error);
        return MSO_STATUS_ERROR;
    }
    return MSO_STATUS_OK;
}

const mso_backend_api kApi = {
    MSO_BACKEND_API_VERSION,
    MSO_BACKEND_KIND_CUDA,
    "cuda",
    &createBackend,
    &destroyBackend,
    &lastError,
    &solveFullGridMin,
};

} // namespace

extern "C" MSO_BACKEND_EXPORT const mso_backend_api* MSO_BACKEND_CALL mso_backend_get_api()
{
    return &kApi;
}
