cmake_minimum_required(VERSION 3.20)

project(MultistageOptimizer LANGUAGES CXX)

option(MSO_BUILD_TESTS "Build unit tests" ON)
option(MSO_BUILD_CUDA_BACKEND "Build CUDA backend plugin" ON)
set(MSO_CUDA_ARCHITECTURES "" CACHE STRING "CUDA architectures for mso_backend_cuda (e.g. native; all-major; 86;90). Defaults to 'native' when supported.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MSO_OUTPUT_DIR "${CMAKE_BINARY_DIR}/bin")

if (WIN32 AND NOT CMAKE_MT)
    file(GLOB MSO_MT_CANDIDATES
        "C:/Program Files (x86)/Windows Kits/10/bin/*/x64/mt.exe"
    )
    list(LENGTH MSO_MT_CANDIDATES MSO_MT_COUNT)
    if (MSO_MT_COUNT GREATER 0)
        list(GET MSO_MT_CANDIDATES 0 MSO_MT_EXE)
        set(CMAKE_MT "${MSO_MT_EXE}" CACHE FILEPATH "Manifest tool" FORCE)
    endif()
endif()

add_executable(MultistageOptimizer
    MultistageOptimizer.cpp
)
target_include_directories(MultistageOptimizer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(MultistageOptimizer PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${MSO_OUTPUT_DIR}")

if (WIN32)
    target_sources(MultistageOptimizer PRIVATE MultistageOptimizer.rc)
endif()

if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    target_compile_options(MultistageOptimizer PRIVATE /W4 /permissive-)
else()
    target_compile_options(MultistageOptimizer PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (MSO_BUILD_CUDA_BACKEND)
    include(CheckLanguage)
    if (WIN32)
        set(_mso_need_msvc_cuda_host_compiler OFF)
        if (NOT CMAKE_CUDA_HOST_COMPILER)
            set(_mso_need_msvc_cuda_host_compiler ON)
        else()
            string(TOLOWER "${CMAKE_CUDA_HOST_COMPILER}" _mso_cuda_host_compiler_lower)
            if (_mso_cuda_host_compiler_lower MATCHES "clang-cl")
                set(_mso_need_msvc_cuda_host_compiler ON)
            endif()
        endif()
    endif()

    if (WIN32 AND _mso_need_msvc_cuda_host_compiler)
        file(GLOB MSO_CL_CANDIDATES
            "C:/Program Files/Microsoft Visual Studio/*/*/VC/Tools/MSVC/*/bin/Hostx64/x64/cl.exe"
            "C:/Program Files (x86)/Microsoft Visual Studio/*/*/VC/Tools/MSVC/*/bin/Hostx64/x64/cl.exe"
        )
        list(LENGTH MSO_CL_CANDIDATES MSO_CL_COUNT)
        if (MSO_CL_COUNT GREATER 0)
            list(GET MSO_CL_CANDIDATES 0 MSO_CL_EXE)
            set(CMAKE_CUDA_HOST_COMPILER "${MSO_CL_EXE}" CACHE FILEPATH "Host compiler for CUDA" FORCE)
        endif()
    endif()
    find_program(MSO_NVCC nvcc HINTS "$ENV{ProgramFiles}/NVIDIA GPU Computing Toolkit/CUDA" PATH_SUFFIXES bin)
    if (MSO_NVCC)
        set(CMAKE_CUDA_COMPILER "${MSO_NVCC}" CACHE FILEPATH "CUDA compiler" FORCE)
    endif()
    check_language(CUDA)
    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        add_library(mso_backend_cuda SHARED
            backends/cuda/mso_backend_cuda.cu
        )
        target_include_directories(mso_backend_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
        set_target_properties(mso_backend_cuda PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED ON
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED ON
            POSITION_INDEPENDENT_CODE ON
            RUNTIME_OUTPUT_DIRECTORY "${MSO_OUTPUT_DIR}"
            LIBRARY_OUTPUT_DIRECTORY "${MSO_OUTPUT_DIR}"
        )
        if (WIN32)
            target_compile_options(mso_backend_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--allow-unsupported-compiler>)
        endif()
        if (MSO_CUDA_ARCHITECTURES STREQUAL "" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
            set(MSO_CUDA_ARCHITECTURES native)
        endif()
        if (NOT MSO_CUDA_ARCHITECTURES STREQUAL "")
            set_property(TARGET mso_backend_cuda PROPERTY CUDA_ARCHITECTURES ${MSO_CUDA_ARCHITECTURES})
        endif()
    else()
        message(STATUS "CUDA compiler not found; skipping mso_backend_cuda.")
    endif()
endif()

if (MSO_BUILD_TESTS)
    enable_testing()

    file(GLOB MSO_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp"
    )
    add_executable(MultistageOptimizerTests ${MSO_TEST_SOURCES})
    target_include_directories(MultistageOptimizerTests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
    )
    set_target_properties(MultistageOptimizerTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${MSO_OUTPUT_DIR}")

    if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        target_compile_options(MultistageOptimizerTests PRIVATE /W4 /permissive-)
    else()
        target_compile_options(MultistageOptimizerTests PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_test(NAME MultistageOptimizerTests COMMAND $<TARGET_FILE:MultistageOptimizerTests>)
endif()
