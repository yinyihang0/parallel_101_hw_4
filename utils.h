#include <chrono>
#include <cmath>
#include <cstdio>

constexpr int N = 48;
constexpr int NumPara = N * 7;

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;

float frand()
{
    return (float)rand() / RAND_MAX * 2 - 1;
}

void printResults(const char* description, double startEnergy, double finalEnergy, long elapsedTime) {
    printf("%-20s Start Energy: %-12.4f Final Energy: %-12.4f Time Elapsed: %ld ms\n",
           description, startEnergy, finalEnergy, elapsedTime);
}

template <class Func>
long benchmark(Func const &func)
{
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}