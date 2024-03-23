#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <iostream>
#include <stdint.h>
#include <string.h>
#include <array>


namespace sse_omp
{
    struct Stars{
        float pxs[N];
        float pys[N];
        float pzs[N];
        float vxs[N];
        float vys[N];
        float vzs[N];
        float masses[N];
    };

    Stars stars;

    void init(const std::array<float, NumPara>& arr){
        for (int i = 0; i < N; i++){
            stars.pxs[i] = arr[i*7];
            stars.pys[i] = arr[i*7+1];
            stars.pzs[i] = arr[i*7+2];
            stars.vxs[i] = arr[i*7+3];
            stars.vys[i] = arr[i*7+4];
            stars.vzs[i] = arr[i*7+5];
            stars.masses[i] = arr[i*7+6];
        }
    }

    void step() {
        size_t len = N;
        float eps2 = eps * eps;
        float gdt = G * dt;
        // #pragma GCC unroll 16
        for(size_t i=0; i<len; i++){
            float dxs[N];
            float dys[N];
            float dzs[N];
            float d2s[N];
            float ivf_d2s[N];
            #pragma omp simd        
            for(size_t j=0; j<len; j++){
                dxs[j] = stars.pxs[j] - stars.pxs[i];
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                dys[j] = stars.pys[j] - stars.pys[i];
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                dzs[j] = stars.pzs[j] - stars.pzs[i];
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                d2s[j] = dxs[j] * dxs[j] + dys[j] * dys[j] + dzs[j] * dzs[j] + eps2;
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                ivf_d2s[j] = 1.0 / (d2s[j] * std::sqrt(d2s[j]));
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                stars.vxs[i] += dxs[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                stars.vys[i] += dys[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
            }
            #pragma omp simd
            for(size_t j=0; j<len; j++){
                stars.vzs[i] += dzs[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
            }
        }

        #pragma omp simd
        for(size_t i=0; i<len; i++){
            stars.pxs[i] += stars.vxs[i] * dt;
        }
        #pragma omp simd
        for(size_t i=0; i<len; i++){
            stars.pys[i] += stars.vys[i] * dt;
        }
        #pragma omp simd
        for(size_t i=0; i<len; i++){
            stars.pzs[i] += stars.vzs[i] * dt;
        }
    }

    float calc() {
        float energy = 0;
        size_t len = N;
        for(size_t i=0; i<len; i++){
            float v2 = stars.vxs[i] * stars.vxs[i] + stars.vys[i] * stars.vys[i] + stars.vzs[i] * stars.vzs[i];
            energy += stars.masses[i] * v2 / 2;
            // #pragma GCC unroll 32
            for(size_t j=0; j<len; j++){
                float dx = stars.pxs[j] - stars.pxs[i];
                float dy = stars.pys[j] - stars.pys[i];
                float dz = stars.pzs[j] - stars.pzs[i];
                float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
                float ivf_d2 = 1.0 / (std::sqrt(d2) * 2);
                energy -= stars.masses[i] * stars.masses[j] * (G * ivf_d2);
            }
        }
        return energy;
    }


};
    

