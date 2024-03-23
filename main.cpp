#include "utils.h"
#include "origin.h"
#include "sse_intrin.h"
#include "sse_omp.h"
#include "sse_intrin_rcp_rsqrt.h"

#include <array>
#include <assert.h>

// initialize global array
std::array<float, NumPara>& init()
{
    // static will only initialize once
    static std::array<float, NumPara> globalArr;

    static bool initialized = false;
    if(!initialized) {
        for (int i = 0; i < NumPara; i++)
        {
            globalArr[i] = frand();
            if(i % 7 == 6) globalArr[i] = 1.0f + frand();
        }
        initialized = true;
    }

    return globalArr;
}


int main(int argc, char const *argv[])
{

    int num = 5000000;

    // get ref, while auto might have to copy
    auto& globalArr = init();

    // init origin and sse
    origin::init(globalArr);
    sse_omp::init(globalArr);
    sse_intrin::init(globalArr);
    sse_intrin_rcp_rsqrt::init(globalArr);


    // 原始版本
    double startEnergy = origin::calc();
    auto dt = benchmark([&]{
        for (int i = 0; i < num; i++) {
            origin::step();
        }
    });
    double finalEnergy = origin::calc();
    printResults("Original:", startEnergy, finalEnergy, dt);

    // SIMD OpenMP 版本
    startEnergy = sse_omp::calc();
    dt = benchmark([&]{
        for (int i = 0; i < num; i++) {
            sse_omp::step();
        }
    });
    finalEnergy = sse_omp::calc();
    printResults("SIMD OpenMP:", startEnergy, finalEnergy, dt);

    // SIMD Intrinsics 版本
    startEnergy = sse_intrin::calc();
    dt = benchmark([&]{
        for (int i = 0; i < num; i++) {
            sse_intrin::step();
        }
    });
    finalEnergy = sse_intrin::calc();
    printResults("SIMD Intrinsics:", startEnergy, finalEnergy, dt);

    // SIMD INTRINSICS DIV 版本
    startEnergy = sse_intrin_rcp_rsqrt::calc();
    dt = benchmark([&]{
        for (int i = 0; i < num; i++) {
            sse_intrin_rcp_rsqrt::step();
        }
    });
    finalEnergy = sse_intrin_rcp_rsqrt::calc();
    printResults("SIMD Intrinsics Div:", startEnergy, finalEnergy, dt);
    return 0;
}
