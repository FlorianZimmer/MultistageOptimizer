cmake_minimum_required(VERSION 3.20)

project(MultistageOptimizer LANGUAGES CXX)

option(MSO_BUILD_TESTS "Build unit tests" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(MultistageOptimizer
    MultistageOptimizer.cpp
)
target_include_directories(MultistageOptimizer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (WIN32)
    target_sources(MultistageOptimizer PRIVATE MultistageOptimizer.rc)
endif()

if (MSVC)
    target_compile_options(MultistageOptimizer PRIVATE /W4 /permissive-)
else()
    target_compile_options(MultistageOptimizer PRIVATE -Wall -Wextra -Wpedantic)
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

    if (MSVC)
        target_compile_options(MultistageOptimizerTests PRIVATE /W4 /permissive-)
    else()
        target_compile_options(MultistageOptimizerTests PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_test(NAME MultistageOptimizerTests COMMAND $<TARGET_FILE:MultistageOptimizerTests>)
endif()
