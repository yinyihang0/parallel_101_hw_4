#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <iostream>
#include <array>
namespace origin
{

    struct Star
    {
        float px, py, pz;
        float vx, vy, vz;
        float mass;
    };

    std::vector<Star> stars;

    void init(const std::array<float, NumPara>& arr)
    {
        for (int i = 0; i < N; i++)
        {
            stars.push_back({arr[i*7], arr[i*7+1], arr[i*7+2], arr[i*7+3], arr[i*7+4], arr[i*7+5], arr[i*7+6]});
        }
    }


    void step()
    {
        for (auto &star : stars)
        {
            for (auto &other : stars)
            {
                float dx = other.px - star.px;
                float dy = other.py - star.py;
                float dz = other.pz - star.pz;
                float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
                d2 *= sqrt(d2);
                star.vx += dx * other.mass * G * dt / d2;
                star.vy += dy * other.mass * G * dt / d2;
                star.vz += dz * other.mass * G * dt / d2;
            }
        }
        for (auto &star : stars)
        {
            star.px += star.vx * dt;
            star.py += star.vy * dt;
            star.pz += star.vz * dt;
        }
    }

    float calc()
    {   
        float energy = 0;
        for (auto &star : stars)
        {
            float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
            energy += star.mass * v2 / 2;
            for (auto &other : stars)
            {
                float dx = other.px - star.px;
                float dy = other.py - star.py;
                float dz = other.pz - star.pz;
                float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
                energy -= other.mass * star.mass * G / sqrt(d2) / 2;
            }
        }
        return energy;
    }

}

