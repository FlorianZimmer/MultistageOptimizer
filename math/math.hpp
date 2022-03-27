#pragma once

#include <iostream>

// Returns value of Binomial Coefficient C(n, k)
// source: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
unsigned long long nCr(int n, int k)
{
    unsigned long long res = 1;

    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;

    // Calculate value of
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (long long i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}

//fill Array with all possible combination aka black magic fuckery
static void createTuple(unsigned char* arr, int sizeX, unsigned long long sizeY, int sum) {
    int max = sum - sizeX + 1;
    arr[0 * sizeY + 0] = max;  //arr[column * sizeY + line]
    for (long i = 1; i < sizeY; i++) {
        if (arr[sizeX * sizeY] == max) {
            return;
        }
        //copy previous line to current
        for (int k = 0; k < sizeX; k++) {
            arr[k * sizeY + i] = arr[k * sizeY + i - 1];
        }
        if (arr[0 * sizeY + i - 1] > 1) {
            arr[0 * sizeY + i] -= 1;
            arr[1 * sizeY + i] += 1;
        }
        else {
            int col = 1;
            while (arr[col * sizeY + i - 1] == 1) {
                col += 1;
            }
            arr[0 * sizeY + i] = arr[col * sizeY + i - 1] - 1;
            arr[(col + 1) * sizeY + i] = arr[(col + 1) * sizeY + i - 1] + 1;
            arr[col * sizeY + i] = 1;
        }

        //print array for line by line. Used for debugging
        //for (int j = 0; j < sizeX; j++) {
        //    int sum = 0;
        //    std::cout << arr[j * sizeY + i] << " ";
        //    sum += arr[j * sizeY + i];
        //}
        //std::cout << sum;
        //std::cout << "\n";
    }
}

//Addition function that detects overflow
unsigned long safeAdd_ULONG(unsigned long lhs, unsigned long rhs)
{
    if (lhs >= 0) {
        if (ULONG_MAX - lhs < rhs) {
            return ULONG_MAX;
        }
    }
    else {
        if (rhs < ULONG_MAX - lhs) {
            return ULONG_MAX;
        }
    }
    return lhs + rhs;
}

unsigned long long safeMULT_ULLONG(unsigned long long a, unsigned long long b) {
    unsigned long long x = a * b;
    if (a != 0 && x / a != b) {
        return ULLONG_MAX;
    }
    else {
        return x;
    }
}

//calculate the maximum precision without causing an overflow when creating the array distributions
int calcMaxPrecUp(int precision, int nStages, unsigned long long maxRAM, unsigned long long nCombinations) {
    int i = precision;
    while (safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCr(i - 1, nStages - 1) * sizeof(unsigned long)) >= maxRAM) {
        unsigned long long temp = safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCr(i - 1, nStages - 1) * sizeof(unsigned long));
        i--;
    }
    return i;
}

//calculate the maximum precision without causing an overflow when creating the array distributions
int calcMaxPrecDown(int precision, int nStages, unsigned long long maxRAM, unsigned long long nCombinations) {
    int i = precision;
    while (safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCombinations * sizeof(unsigned long)) >= maxRAM) {
        unsigned long long temp = safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCombinations * sizeof(unsigned long));
        i--;
    }
    i--;
    return i;
}

static double* ArrToPerc(double* newArr, unsigned char* arr, int sizeX, int precision) {
    for (int i = 0; i < sizeX; i++) {
        newArr[i] = ((double)arr[i] / precision) * 100;
    }
    return newArr;
}