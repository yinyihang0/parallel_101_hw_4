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


namespace sse_intrin_rcp_rsqrt
{

    void print128_num(__m128 var)
    {
        float val[4];
        memcpy(val, &var, sizeof(val));
        printf("Numerical: %f %f %f %f \n",
               val[0], val[1], val[2], val[3]);
    }

    struct SIMDStar
    {
        alignas(16) __m128 px, py, pz;
        alignas(16) __m128 vx, vy, vz;
        alignas(16) __m128 mass;
    };


    // structure of array
    template <int NUM>
    class Stars
    {
    public:
        alignas(16) float px[NUM], py[NUM], pz[NUM];
        alignas(16) float vx[NUM], vy[NUM], vz[NUM];
        alignas(16) float mass[NUM];

        SIMDStar extractSIMDStar(int index)
        {
            // static_assert(index % 4 == 0, "index must be multiple of 4");
            return {
                _mm_load_ps(px + index), _mm_load_ps(py + index), _mm_load_ps(pz + index),
                _mm_load_ps(vx + index), _mm_load_ps(vy + index), _mm_load_ps(vz + index),
                _mm_load_ps(mass + index)
            };
        }
        SIMDStar replicateStar(int index)
        {
            return {
                _mm_set_ps1(px[index]), _mm_set_ps1(py[index]), _mm_set_ps1(pz[index]),
                _mm_set_ps1(vx[index]), _mm_set_ps1(vy[index]), _mm_set_ps1(vz[index]),
                _mm_set_ps1(mass[index])
            };
        }
    };

    Stars<N> stars;

    void init(const std::array<float, NumPara>& arr)
    {
        for (int i = 0; i < 48; i++)
        {
            stars.px[i] = arr[i*7];
            stars.py[i] = arr[i*7+1];
            stars.pz[i] = arr[i*7+2];
            stars.vx[i] = arr[i*7+3];
            stars.vy[i] = arr[i*7+4];
            stars.vz[i] = arr[i*7+5];
            stars.mass[i] = arr[i*7+6];
        }
    }

    // sum all four elements in __m128
    float sum_m128(__m128 var)
    {
        __m128 var1 = _mm_hadd_ps(var, var);
        __m128 var2 = _mm_hadd_ps(var1, var1);
        return _mm_cvtss_f32(var2);
    }

    __m128 eps_squares = _mm_set_ps1(eps * eps);
    __m128 ones = _mm_set_ps1(1.0f);
    __m128 G_times_dt = _mm_set_ps1(G * dt);
    __m128 dts = _mm_set_ps1(dt);
    __m128 zeros = _mm_set_ps1(0.0f);
    

    int num = 0;
    void step()
    {
        for (int i = 0; i < 48; i++)
        {
            auto star = stars.replicateStar(i);
            __m128 vx_update = zeros;
            __m128 vy_update = zeros;
            __m128 vz_update = zeros;

            for (int j = 0; j < 48; j += 4)
            {
                auto other = stars.extractSIMDStar(j);
                __m128 dx = _mm_sub_ps(other.px, star.px);
                __m128 dy = _mm_sub_ps(other.py, star.py);
                __m128 dz = _mm_sub_ps(other.pz, star.pz);
            

                __m128 d2 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz)), eps_squares);

                __m128 inv_d2 = _mm_rcp_ps(d2);
                // eliminate the denpendency: not denpend on inv_d2
                __m128 inv_d = _mm_rsqrt_ps(d2);
                __m128 d3_inv = _mm_mul_ps(inv_d2, inv_d);


                vx_update = _mm_add_ps(vx_update, _mm_mul_ps(_mm_mul_ps(dx, other.mass), _mm_mul_ps(d3_inv, G_times_dt)));
                vy_update = _mm_add_ps(vy_update, _mm_mul_ps(_mm_mul_ps(dy, other.mass), _mm_mul_ps(d3_inv, G_times_dt)));
                vz_update = _mm_add_ps(vz_update, _mm_mul_ps(_mm_mul_ps(dz, other.mass), _mm_mul_ps(d3_inv, G_times_dt)));
            }
            
            *(stars.vx + i) += sum_m128(vx_update);
            *(stars.vy + i) += sum_m128(vy_update);
            *(stars.vz + i) += sum_m128(vz_update);
        }
      
        for (int i = 0; i < 48; i += 4)
        {
            auto star = stars.extractSIMDStar(i);
            _mm_store_ps(stars.px + i, _mm_add_ps(star.px, _mm_mul_ps(star.vx, dts)));
            _mm_store_ps(stars.py + i, _mm_add_ps(star.py, _mm_mul_ps(star.vy, dts)));
            _mm_store_ps(stars.pz + i, _mm_add_ps(star.pz, _mm_mul_ps(star.vz, dts)));
        }
    }


    __m128 halfs = _mm_set_ps1(0.5f);
    __m128 doublehalfs = _mm_set_ps1(0.25f);
    __m128 Gs = _mm_set_ps1(G);
    float calc()
    {
        __m128 energy = _mm_set_ps1(0.0f);
        for (int i = 0; i < 48; i++)
        {
            auto star = stars.replicateStar(i);
            __m128 v2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(star.vx, star.vx), _mm_mul_ps(star.vy, star.vy)), _mm_mul_ps(star.vz, star.vz));

            energy = _mm_add_ps(energy, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(star.mass, v2), halfs), doublehalfs));

            for (int j = 0; j < 48; j += 4)
            {
                auto other = stars.extractSIMDStar(j);
                __m128 dx = _mm_sub_ps(other.px, star.px);
                __m128 dy = _mm_sub_ps(other.py, star.py);
                __m128 dz = _mm_sub_ps(other.pz, star.pz);
                __m128 d2 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz)), eps_squares);
                __m128 sqrt_d2_2 = _mm_div_ps(ones, _mm_sqrt_ps(d2));
                energy = _mm_sub_ps(energy, _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(other.mass, star.mass), sqrt_d2_2), Gs), halfs));
            }
        }
        return (sum_m128(energy));
    }
}
