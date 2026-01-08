#pragma once

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#define MSO_BACKEND_CALL __cdecl
#define MSO_BACKEND_EXPORT __declspec(dllexport)
#else
#define MSO_BACKEND_CALL
#define MSO_BACKEND_EXPORT __attribute__((visibility("default")))
#endif

#define MSO_BACKEND_API_VERSION 1u

typedef enum mso_status {
    MSO_STATUS_OK = 0,
    MSO_STATUS_ERROR = 1,
    MSO_STATUS_INVALID_ARGUMENT = 2,
    MSO_STATUS_UNSUPPORTED = 3,
    MSO_STATUS_NO_DEVICE = 4,
    MSO_STATUS_OUT_OF_MEMORY = 5,
} mso_status;

typedef enum mso_backend_kind {
    MSO_BACKEND_KIND_UNKNOWN = 0,
    MSO_BACKEND_KIND_CUDA = 1,
    MSO_BACKEND_KIND_ROCM = 2,
} mso_backend_kind;

typedef struct mso_backend_create_params {
    int32_t device_id; // -1 = default device
    uint32_t flags;    // reserved, must be 0
} mso_backend_create_params;

typedef struct mso_full_grid_problem {
    uint32_t stage_count;
    uint32_t precision;
    uint64_t combinations; // expected combinations count (may be 0 to let backend compute)
    double payload;
    const double* inert_mass_rest;        // [stage_count]
    const double* inert_mass_fract;       // [stage_count]
    const double* one_minus_inert_fract;  // [stage_count]
    const double* inert_coeff;            // [stage_count]
    const double* ratio_tables_flat;      // [stage_count * (precision + 1)], stage-major
} mso_full_grid_problem;

typedef struct mso_full_grid_solution {
    double best_mass;
    uint64_t best_index;
    uint8_t* best_units;     // optional output buffer [best_units_len]
    uint32_t best_units_len; // must be >= stage_count when best_units != NULL
} mso_full_grid_solution;

typedef struct mso_backend_api {
    uint32_t api_version;
    uint32_t backend_kind;
    const char* backend_name;

    mso_status(MSO_BACKEND_CALL* create)(const mso_backend_create_params* params, void** out_backend);
    void(MSO_BACKEND_CALL* destroy)(void* backend);
    const char*(MSO_BACKEND_CALL* last_error)(void* backend);

    mso_status(MSO_BACKEND_CALL* solve_full_grid_min)(void* backend,
                                                      const mso_full_grid_problem* problem,
                                                      mso_full_grid_solution* solution);
} mso_backend_api;

typedef const mso_backend_api*(MSO_BACKEND_CALL* mso_backend_get_api_fn)();

