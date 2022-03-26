/*
Filename: Multistage_optimizer_0.1_CPP.cpp
Description : Calculate optimal mass for each stage given a DeltaV requirement on a multistage Rocket
Source of math equations : http://www.projectrho.com/public_html/rocket/multistage.php
Author: Florian Zimmer
Version : 0.4
Last Changed : 17.03.2021
Changelog :
0.1 basic functionality/port from python project
0.2 debugging
0.3 json input
0.4 cleanup code
0.5 fix Perfomance Issues and Multithreading
*/

#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include "json.hpp"
using json = nlohmann::json;

using namespace std;

class Engine {
public:
    string name;
    double isp_sl;      
    double isp_vac;
    double mass;        //kg
    Engine(string name, double isp_sl, double isp_vac, double mass) {
        this->name = name;
        this->isp_sl = isp_sl;
        this->isp_vac = isp_vac;
        this->mass = mass;
    }
    Engine() {
        this->name = "null";
        this->isp_sl = 0;
        this->isp_vac = 0;
        this->mass = 0;
    }
};

class Stage {
public:
    int num_engines;
    Engine engine;
    double tot_mass;
    double wet_mass;
    double dry_mass;
    double prop_mass;  //calc
    double inert_mass_tankonly; //calc
    double inert_mass_rest; //calc
    double inert_mass_fract; //calc
    double payload;
    Stage(double payload, Engine engine, int num_engines, double tot_mass, double dry_mass, double wet_mass) {
        this->payload = payload;

        //can be deleted as well as not necessary info. Weight can be calculated as seen in inert_mass_rest
        this->num_engines = num_engines;
        this->engine = engine;
        this->tot_mass = tot_mass;
        this->wet_mass = wet_mass;
        this->dry_mass = dry_mass;
        prop_mass = wet_mass - dry_mass;
        //calculations are more precise, by knowing what is a scalable inert mass and what not, when increasing or decreasing the fuel in a tank, 
        //that is why the to different variables for inert masses exist, one for calculation of inert mass fraction and the other non scalable, which can include the mass of engine(s), rcs, control unit, etc.
        inert_mass_rest = tot_mass - wet_mass;

        //can be deleted -> dry mass
        inert_mass_tankonly = dry_mass;     //tot_mass - prop_mass - inert_mass_rest;   //inert mass of the fuel tank only
        inert_mass_fract = (double)inert_mass_tankonly / wet_mass;      
    };                                                                  
};

class Rocket {
public:
    vector<Stage> stages;
    int deltaV;
    double g_body;
    Rocket(double g_body, double mission_payload, int deltaV, vector<Stage> stages) //: stages(nStages)
    {
        this->g_body = g_body;
        this->deltaV = deltaV;
        for (int i = 0; i < stages.size(); i++) {
            if (i == 0) {
                Stage newStage = Stage(mission_payload, stages[i].engine, stages[i].num_engines, stages[i].tot_mass, stages[i].dry_mass, stages[i].wet_mass);
                this->stages.push_back(newStage);
            }
            else
            {
                Stage newStage = Stage(this->stages[i - 1].tot_mass  + this->stages[i - 1].payload, stages[i].engine, stages[i].num_engines, stages[i].tot_mass, stages[i].dry_mass, stages[i].wet_mass);
                this->stages.push_back(newStage);
            }
        }
    }
};

static json readJson(string path) {
    std::ifstream ifs(path);
    json j = json::parse(ifs);
    return j;
}

static vector<Engine> find_engines_on_Stage(vector<Engine> engineList, json file) {
    vector<Engine> engines_on_Stage;
    for (int i = 0; i < file["stages"].size(); i++)
    {
        for (int j = 0; j < engineList.size(); j++) {
            if (file["stages"][i]["engine"] == engineList[j].name) {
                engines_on_Stage.push_back(engineList[j]);
            }
        }
    }
    return engines_on_Stage;
}

static vector<Stage> extractStages(vector<Engine> engineList, json file) {
    vector<Stage> all_stages;
    vector<Engine> engines = find_engines_on_Stage(engineList, file);
    for (int i = 0; i < file["stages"].size(); i++) {
        Stage tempStage(0, engines[i], file["stages"][i]["engine_count"], file["stages"][i]["total_mass"], file["stages"][i]["dry_mass"], file["stages"][i]["wet_mass"]);
        all_stages.push_back(tempStage);
    }
    return all_stages;
}

static vector<Engine> importEngines(string path) {
    vector<Engine> engineList;
    json file = readJson(path);
    for (int i = 0; i < file["engines"].size(); i++) {
        engineList.push_back(Engine(file["engines"][i]["name"], file["engines"][i]["isp_sl"], file["engines"][i]["isp_vac"], file["engines"][i]["mass"]));
    }
    return engineList;
}

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

static double* ArrToPerc(double* newArr, unsigned char* arr, int sizeX, int precision) {
    for (int i = 0; i < sizeX; i++) {
        newArr[i] = ((double)arr[i]/precision)*100;
    }
    return newArr;
}

//calculate factorial of number n
static double fact(int n)
{
    if (n <= 1)
    {
        return n;
    }
    else { return n * fact(n - 1); }
}

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

//mass ratio
static double massRatio(double deltaVfraction, int deltaVmission, double gravity, double isp) {
    return std::exp((deltaVfraction * deltaVmission) / (gravity * isp));
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


static unsigned long calcMass(unsigned long payloadMass, float inertMassFract, float deltaVfraction, int deltaVmission, float gravity, float isp, double inert_mass_rest) {
    double R = massRatio(deltaVfraction, deltaVmission, gravity, isp);
    unsigned long mass = 0;

    //when inertMassFract * R >= 1 than the targeted deltav is phisically not achievable with the given effiency of the engine and the inert mass fraction of the tank
    if (inertMassFract * R >= 1){
        mass = ULONG_MAX; 
    }
    else {
        //safe addition as overflows occure here regularly, if overflow is detected the result is ULONG_MAX
        double propMass = ( (safeAdd_ULONG(payloadMass, inert_mass_rest) * ((R - 1) * (1 - inertMassFract)) ) ) / (1 - (inertMassFract * R));   //adding non-scalable inert mass as payload
        unsigned long inertMass = safeAdd_ULONG(((inertMassFract / (1 - inertMassFract)) * propMass), inert_mass_rest);     
        mass = safeAdd_ULONG(safeAdd_ULONG(propMass, inertMass), payloadMass);
    }

    return mass;
}

//calculate total mass of rocket for all distributions
static void distMass(bool verbose, int start, int end, long long nCombinations, Rocket rocket, unsigned long * massDistributions, unsigned char* distributions, int precision) {
    int min = INT_MAX;

    //getting all constant vars of the rocket here as to not repeatedly call them over an over for every iteration in for-loop -> speed ups program by about 30%
    int rocket_size = rocket.stages.size();
    int deltaV = rocket.deltaV;
    double g_body = rocket.g_body;
    double* inert_mass_rest = new double[rocket_size];
    double* inert_mass_fract = new double[rocket_size];
    double* isp_vac = new double[rocket_size];
    for (int j = 0; j < rocket_size; j++) {
        inert_mass_rest[j] = rocket.stages[j].inert_mass_rest;
        inert_mass_fract[j] = rocket.stages[j].inert_mass_fract;
        isp_vac[j] = rocket.stages[j].engine.isp_vac;
    }

    for (int i = start; i < end; i++) {
        //if (i == 258100748){
        //        __debugbreak();
        //}
        for (int j = 0; j < rocket_size; j++) {
            if (j == 0) {
                massDistributions[i] = calcMass(rocket.stages[j].payload, inert_mass_fract[j], distributions[j * nCombinations + i] / (double)precision, deltaV, g_body, isp_vac[j], inert_mass_rest[j]);    //uppermost stage only has mission payload as payload
            }
            //else if (j == rocket.stages.size() - 1) {
            //    massDistributions[i] += calcMass(rocket.stages[j].payload, rocket.stages[j].inert_mass_fract, distributions[j * nCombinations + i] / (double)precision, rocket.deltaV, rocket.g_body, rocket.stages[j].engine.isp_sl, rocket.stages[j]);    //first stage assumes surface level isp
            //}
            else {
                massDistributions[i] = calcMass(massDistributions[i], inert_mass_fract[j], distributions[j * nCombinations + i] / (double)precision, deltaV, g_body, isp_vac[j], inert_mass_rest[j]);
            }
        }
        //output for debbugging
        //only outputs into console, when new best mass is found
        //only coherent with single-threaded execution
        if (verbose == true) {
            if (massDistributions[i] < min) {
                min = massDistributions[i];
                std::cout << "Distribution number: " << i << "\t" << massDistributions[i] << "\n";
                for (int j = 0; j < rocket.stages.size(); j++) {
                    std::cout << "Stage " << rocket_size - j << ": " << (int)distributions[j * nCombinations + i] << "\t";
                }
                std::cout << "\n\n";
            }
        }
    }
    std::cout << "thread " << (int)((double)((double)end / (double)nCombinations)*24) << " closed\n";
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
//calculate the maximum precision without causing an overflow when creating the array distributions
int calcMaxPrecUp(int precision, int nStages, unsigned long long maxRAM, unsigned long long nCombinations) {
    int i = precision;
    while (safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCr(i - 1, nStages - 1) * sizeof(unsigned long)) >= maxRAM) {
        unsigned long long temp = safeMULT_ULLONG(nCr(i - 1, nStages - 1), nStages) + (nCr(i - 1, nStages - 1) * sizeof(unsigned long));
        i--;
    }
    return i;
}

int main()
{
    //import Data from json files
    vector<Engine> engineList = importEngines("engines.json");
    Rocket rocket = importRocket("rocket.json", engineList);

    //switches -> should be moved to function call or config file
    bool verbose = false;       //output for debbugging
    bool useMultiCore = true;  //using multicore and verbose at the same time can cause gibberish output
    //bool highPrec = true;
    unsigned long long maxRAM = pow(2, 34); //max allowed RAM usage, is the limiting size for the array distributions -> 2^34 = 16GiB

    //max is 255 as values are stored as 1 byte chars in the array distributions
    int precision = 72;        //increment of distributions for calculations; too high numbers can cause overflow; should warn user of that -> in the future

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