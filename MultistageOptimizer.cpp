/*
Filename: MultistageOptimizer.cpp
Description : Calculate optimal mass for each stage given a DeltaV requirement on a multistage Rocket
Source of math equations : http://www.projectrho.com/public_html/rocket/multistage.php
Author: Florian Zimmer
Version : 0.6
Last Changed : 26.03.2021
Changelog :
0.1 basic functionality/port from python project
0.2 debugging
0.3 json input
0.4 cleanup code
0.5 fix Perfomance Issues and Multithreading
0.6 split files
0.7 add config
*/

#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include "read_json/json.hpp"
using json = nlohmann::json;

#include "rocket_definition/Stage.hpp"
#include "rocket_definition/Engine.hpp"
#include "rocket_definition/Rocket.hpp"
#include "math/math.hpp"
#include "math/massCalc.hpp"
#include "read_json/read_json.hpp"

using namespace std;

static Rocket importRocket(string path, vector<Engine> engineList) {
    json file = readJson(path);
    double g_body = file["g_body"];
    Rocket rocket(file["g_body"], file["mission_payload"], file["deltaVmission"], extractStages(engineList, file));
    return rocket;
}

static void output2Darr(unsigned char* arr, int sizeY, int sizeX)
{
    for (int i = 0; i < sizeY; i++) {
        for (int j = 0; j < sizeX; j++) {
            std::cout << (int)arr[j * sizeY + i] << "\t";
        }
        std::cout << "\n";
    }
}

static void outputArr(double* arr, int sizeX)
{
    for (int i = 0; i < sizeX; i++) {
        std::cout << arr[i] << "%\t";
    }
}

int main()
{
    //import Data from json files
    vector<Engine> engineList = importEngines("config/engines.json");
    Rocket rocket = importRocket("config/rocket.json", engineList);

    //switches -> should be moved to function call or config file
    bool verbose = false;       //output for debbugging
    bool useMultiCore = true;  //using multicore and verbose at the same time can cause gibberish output
    //bool highPrec = true;
    unsigned long long maxRAM = pow(2, 34); //max allowed RAM usage, is the limiting size for the array distributions -> 2^34 = 16GiB

    //max is 255 as values are stored as 1 byte chars in the array distributions
    int precision = 58;        //increment of distributions for calculations; too high numbers can cause overflow; should warn user of that -> in the future

    unsigned long long nCombinations = nCr(precision - 1, rocket.stages.size() - 1); //stars and bars theorem one

    //check if the calculations will stay within the max allowed RAM usage
    if (safeMULT_ULLONG(nCombinations, rocket.stages.size()) + (nCombinations * sizeof(unsigned long)) >= maxRAM) {
        std::cout << "Precision too high for the amount of stages!\n\nDecrease precision, number of stages or increase allowed RAM until this error dissappears. " << "\n\n";
        std::cout << "For the number of stages and RAM currently given, the maximum precision can be: " << calcMaxPrecUp(precision, rocket.stages.size(), maxRAM, nCombinations) << "\n\n";
        return -1;
    }
    else {
        std::cout << "For the number of stages currently given, the maximum precision can be raised to: " << calcMaxPrecDown(255, rocket.stages.size(), maxRAM, nCombinations) << "\n\n";

        std::cout << "Combinations: " << nCombinations << "\n\n";

        //create static 2D array using unsigned char, to save memory, because bigger numbers than 255 (even far lower numbers for that matter) for precision are unnecessarily computaionally expensive
        //size of array cant be bigger than (2^32) - 1  4294967296
        //unsigned char* test = new unsigned char[pow(2,36)];
        unsigned char* distributions = new unsigned char[nCombinations * rocket.stages.size()]; //access with distributions[j*nCombinations+i]

        //timer for function runtime checking
        std::clock_t start;
        double duration;
        start = std::clock();

        //fill Array with 1s for application of stars and bars theorem one
        std::fill_n(distributions, nCombinations * rocket.stages.size(), 1);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << duration << " seconds to Fill Array with 1s" << '\n';

        //timer for function runtime checking
        start = std::clock();

        createTuple(distributions, rocket.stages.size(), nCombinations, precision);

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << duration << " seconds to create the tuple of all possible distributions\n\n";

        //print all possible fractions --> Debugging
        //output2Darr(distributions, nCombinations, rocket.stages.size());

        //create array with double formatted (/100) fractions so 1 = 0.01 and 50 = 0.5

        //create array for results
        unsigned long* massDistributions = new unsigned long[nCombinations];
        //fill Array with 0s
        std::fill_n(massDistributions, nCombinations, 0);

        //timer for function runtime checking
        start = std::clock();

        //find out how many threads to create
        int nThreads;
        if (useMultiCore == true) {
            cout << "using " << std::thread::hardware_concurrency() << " threads\n";
            nThreads = std::thread::hardware_concurrency();
        }
        else
        {
            cout << "using 1 thread\n";
            nThreads = 1;
        }

        long long threadBlock = (nCombinations / nThreads);
        std::vector<std::thread> threads;

        //call function to calculate all possible distributions of masses in nTreads number of threads
        for (int k = 0, begin = 0, end = threadBlock; k < nThreads; k++, begin += threadBlock, end += threadBlock) {
            threads.push_back(
                std::thread([verbose, begin, end, nCombinations, rocket, massDistributions, distributions, precision] {
                    distMass(verbose, begin, end, nCombinations, rocket, massDistributions, distributions, precision);
                    })
            );
        }

        //wait for all threads to complete their tasks
        for (auto& th : threads) {
            th.join();
        }

        std::cout << "calculating rest" << "\n\n";

        //in case of a number of possibilitys not evenly divisible by number of threads do the rest of them single threaded here
        for (int i = threadBlock * nThreads; i < nCombinations; i++) {
            distMass(verbose, i, nCombinations, nCombinations, rocket, massDistributions, distributions, precision);
        }

        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << duration << " seconds to calculate mass off all possible distributions" << "\n\n";

        //timer for function runtime checking
        start = std::clock();

        //find smallest calculated mass
        int min = INT_MAX;
        long long index = 0;
        unsigned char* bestDistro = new unsigned char[rocket.stages.size()];
        for (int i = 0; i < nCombinations; i++) {
            if (massDistributions[i] < min) {
                min = massDistributions[i];
                index = i;
                for (int j = 0; j < rocket.stages.size(); j++) {
                    int tmp = distributions[j * nCombinations + i];
                    bestDistro[j] = distributions[j * nCombinations + i];
                }
            }
        }


        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << duration << " seconds to find minimum mass" << "\n";

        std::cout << min << " kg at Distribution number " << index << "\n";
        std::cout << (min / 1000) << " t" << "\n\n";
        std::cout << "best deltaV distribution is:\n";
        double* newArr = new double[rocket.stages.size()];
        newArr = ArrToPerc(newArr, bestDistro, rocket.stages.size(), precision);
        outputArr(newArr, rocket.stages.size());
        std::cout << "\n";
        for (int i = 0; i < rocket.stages.size(); i++)
        {
            std::cout << "Stage " << rocket.stages.size() - i << ":\t" << ((((double)bestDistro[i]) / precision) * rocket.deltaV) << "\n";  //delatV of each stage
        }


        delete[] distributions;
    }
    
}