# Intrinsic Optimization

[中文版](./README.zh.md)

## Project Introduction

Optimization of celestial system energy calculation methods using processor built-in functions (intrinsics). Hotspot analysis and microarchitecture analysis were performed using Intel VTune Amplifier, further improving the algorithm's execution efficiency.

## Environment Requirements

- C++ compiler (supporting C++11 standard)
- Linux environment toolchain (recommended to use WSL with VsCode for programming)

## Specific Implementation

- **SSE and OpenMP version (`sse_omp.h`)**: By separating the calculation of `x`, `y`, `z` and utilizing OpenMP for parallel computing, the program's execution efficiency was improved.

- **Fully using Intrinsic version (`sse_intrin.h`)**: Going deeper, directly calling processor built-in functions, bypassing the performance overhead of high-level languages, further optimized data processing.

- **Optimized with fast reciprocal and square root instructions (`sse_intrin_rcp_rsqrt.h`)**: After conducting Intel VTune analysis targeting performance bottlenecks, we found that division operations became the main limiting factor, leading to a decrease in instruction execution efficiency (most instructions were waiting for the result of division). To address this issue, we introduced fast reciprocal (rcp) and fast square root (rsqrt) instructions to replace traditional division and square root operations. This improvement significantly enhanced the program's performance, effectively reducing execution blocking caused by division operations, thus making the instruction execution smoother and improving the overall running efficiency of the program.

## Results
Results from the analysis using Intel VTune:
Microarchitecture Exploration:
![Alt text](image.png)

HPC Performance Characterization
![Alt text](image-1.png)
