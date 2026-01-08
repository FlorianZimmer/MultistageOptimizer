#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace testfw {

struct Failure final : public std::exception {
    std::string message;

    explicit Failure(std::string message)
        : message(std::move(message))
    {
    }

    const char* what() const noexcept override
    {
        return message.c_str();
    }
};

using TestFn = void (*)();

struct TestCase final {
    const char* name = "";
    TestFn fn = nullptr;
};

inline std::vector<TestCase>& registry()
{
    static std::vector<TestCase> tests;
    return tests;
}

struct Registrar final {
    Registrar(const char* name, TestFn fn)
    {
        registry().push_back(TestCase{name, fn});
    }
};

inline std::string formatLocation(const char* file, int line)
{
    std::ostringstream oss;
    oss << file << ":" << line;
    return oss.str();
}

inline void require(bool condition, const char* expression, const char* file, int line)
{
    if (condition) {
        return;
    }
    std::ostringstream oss;
    oss << formatLocation(file, line) << ": REQUIRE failed: " << expression;
    throw Failure(oss.str());
}

inline void requireMsg(bool condition, const char* expression, std::string_view message, const char* file, int line)
{
    if (condition) {
        return;
    }
    std::ostringstream oss;
    oss << formatLocation(file, line) << ": REQUIRE failed: " << expression << " (" << message << ")";
    throw Failure(oss.str());
}

inline bool stringContains(std::string_view haystack, std::string_view needle)
{
    return haystack.find(needle) != std::string_view::npos;
}

inline bool near(double actual, double expected, double epsilon)
{
    if (!std::isfinite(actual) || !std::isfinite(expected)) {
        return false;
    }
    double diff = std::fabs(actual - expected);
    double scale = std::max(1.0, std::max(std::fabs(actual), std::fabs(expected)));
    return diff <= epsilon * scale;
}

inline int runAll()
{
    int failed = 0;
    const auto& tests = registry();

    for (const auto& test : tests) {
        if (!test.fn) {
            continue;
        }
        try {
            test.fn();
        }
        catch (const Failure& e) {
            failed++;
            std::cerr << "[FAIL] " << test.name << "\n" << e.what() << "\n";
        }
        catch (const std::exception& e) {
            failed++;
            std::cerr << "[FAIL] " << test.name << "\nUnhandled std::exception: " << e.what() << "\n";
        }
        catch (...) {
            failed++;
            std::cerr << "[FAIL] " << test.name << "\nUnhandled non-std exception\n";
        }
    }

    if (failed == 0) {
        std::cout << "[PASS] " << tests.size() << " test(s)\n";
        return 0;
    }
    std::cerr << "[FAIL] " << failed << " of " << tests.size() << " test(s)\n";
    return 1;
}

} // namespace testfw

#define TESTFW_CONCAT_INNER(a, b) a##b
#define TESTFW_CONCAT(a, b) TESTFW_CONCAT_INNER(a, b)

#define TEST_CASE(name) TESTFW_TEST_CASE_IMPL(name, __COUNTER__)
#define TESTFW_TEST_CASE_IMPL(name, id)         \
    static void TESTFW_CONCAT(test_, id)();     \
    static testfw::Registrar TESTFW_CONCAT(     \
        test_registrar_,                        \
        id)(name, &TESTFW_CONCAT(test_, id));   \
    static void TESTFW_CONCAT(test_, id)()

#define REQUIRE(expr) testfw::require((expr), #expr, __FILE__, __LINE__)
#define REQUIRE_MSG(expr, message) testfw::requireMsg((expr), #expr, (message), __FILE__, __LINE__)

#define CHECK(expr) REQUIRE(expr)
#define CHECK_MSG(expr, message) REQUIRE_MSG(expr, message)

#define CHECK_NEAR(actual, expected, epsilon)                                                   \
    do {                                                                                        \
        const double testfw_actual = static_cast<double>(actual);                               \
        const double testfw_expected = static_cast<double>(expected);                           \
        REQUIRE_MSG(testfw::near(testfw_actual, testfw_expected, static_cast<double>(epsilon)), \
                    "actual=" + std::to_string(testfw_actual) +                                 \
                        " expected=" + std::to_string(testfw_expected) +                        \
                        " eps=" + std::to_string(static_cast<double>(epsilon)));                \
    } while (0)

#define CHECK_THROWS(expr)                                    \
    do {                                                      \
        bool testfw_threw = false;                            \
        try {                                                 \
            (void)(expr);                                     \
        }                                                     \
        catch (...) {                                         \
            testfw_threw = true;                              \
        }                                                     \
        REQUIRE_MSG(testfw_threw, "Expected exception: " #expr); \
    } while (0)

#define CHECK_THROWS_WITH(expr, substring)                            \
    do {                                                              \
        bool testfw_threw = false;                                    \
        std::string testfw_message;                                   \
        try {                                                         \
            (void)(expr);                                             \
        }                                                             \
        catch (const std::exception& e) {                             \
            testfw_threw = true;                                      \
            testfw_message = e.what();                                \
        }                                                             \
        catch (...) {                                                 \
            testfw_threw = true;                                      \
        }                                                             \
        REQUIRE_MSG(testfw_threw, "Expected exception: " #expr);      \
        REQUIRE_MSG(testfw::stringContains(testfw_message, substring), \
                    "Expected message to contain '" + std::string(substring) + "', got '" + testfw_message + "'"); \
    } while (0)
