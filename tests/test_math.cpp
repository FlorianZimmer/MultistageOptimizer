#include <climits>
#include <cstddef>
#include <limits>
#include <vector>

#include "test_framework.hpp"

#include "math/math.hpp"

TEST_CASE("nCr: handles invalid inputs")
{
    REQUIRE(nCr(-1, 0) == 0);
    REQUIRE(nCr(5, -1) == 0);
    REQUIRE(nCr(5, 6) == 0);
}

TEST_CASE("nCr: basic values and symmetry")
{
    REQUIRE(nCr(0, 0) == 1);
    REQUIRE(nCr(5, 0) == 1);
    REQUIRE(nCr(5, 5) == 1);
    REQUIRE(nCr(5, 2) == 10);
    REQUIRE(nCr(5, 3) == 10);

    for (int n = 0; n <= 20; n++) {
        for (int k = 0; k <= n; k++) {
            REQUIRE(nCr(n, k) == nCr(n, n - k));
        }
    }
}

TEST_CASE("nCr: clamps large values to ULLONG_MAX on overflow")
{
    REQUIRE(nCr(100, 50) == ULLONG_MAX);
}

TEST_CASE("safeAdd_ULLONG: saturates on overflow")
{
    REQUIRE(safeAdd_ULLONG(1, 2) == 3);
    REQUIRE(safeAdd_ULLONG(ULLONG_MAX, 1) == ULLONG_MAX);
    REQUIRE(safeAdd_ULLONG(ULLONG_MAX - 1, 2) == ULLONG_MAX);
}

TEST_CASE("safeMULT_ULLONG: saturates on overflow")
{
    REQUIRE(safeMULT_ULLONG(2, 3) == 6);
    REQUIRE(safeMULT_ULLONG(ULLONG_MAX, 2) == ULLONG_MAX);
    REQUIRE(safeMULT_ULLONG(0, ULLONG_MAX) == 0);
}

TEST_CASE("requiredBytesForPrecision: matches manual sizing for small values")
{
    const int precision = 10;
    const int stages = 2;
    const std::size_t massTypeSize = sizeof(double);

    const unsigned long long combos = nCr(precision - 1, stages - 1); // nCr(9, 1) = 9
    REQUIRE(combos == 9);

    const unsigned long long expected =
        static_cast<unsigned long long>(combos * stages) +
        static_cast<unsigned long long>(combos * massTypeSize);

    REQUIRE(requiredBytesForPrecision(precision, stages, massTypeSize) == expected);
}

TEST_CASE("calcMaxPrecUp: decreases precision until required bytes fit")
{
    const int precision = 10;
    const int stages = 2;
    const unsigned long long maxRAM = 90; // requiredBytesForPrecision(10,2,8) == 90
    const std::size_t massTypeSize = sizeof(double);

    const int max_prec = calcMaxPrecUp(precision, stages, maxRAM, massTypeSize);
    REQUIRE(max_prec == 9);
    REQUIRE(requiredBytesForPrecision(max_prec, stages, massTypeSize) < maxRAM);
}

TEST_CASE("ArrToPerc: converts distribution units to percentages")
{
    unsigned char arr[] = { 1, 2, 7 };
    double out[] = { 0.0, 0.0, 0.0 };

    ArrToPerc(out, arr, 3, 10);
    CHECK_NEAR(out[0], 10.0, 1e-12);
    CHECK_NEAR(out[1], 20.0, 1e-12);
    CHECK_NEAR(out[2], 70.0, 1e-12);
}

