/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common.hpp!"
#endif

#include <thread>

namespace bhc {

// #define BHC_USE_HIGH_PRIORITY_THREADS 1

void SetupThread();

inline uint32_t GetNumCores(bool &singlethread)
{
#ifdef BHC_BUILD_CUDA
    if(singlethread) {
        GlobalLog("Single threaded mode is nonsense on CUDA, ignoring");
        singlethread = false;
    }
#endif
    uint32_t cores = std::thread::hardware_concurrency();
    if(singlethread || cores < 1u) cores = 1u;
#ifdef BHC_USE_HIGH_PRIORITY_THREADS
    if(cores > 1) --cores; // Leave one core unused to avoid locking up system
#endif
    return cores;
}

class Stopwatch {
public:
    Stopwatch() {}
    inline void tick() { tstart = std::chrono::high_resolution_clock::now(); }
    inline void tock(const char *label)
    {
        using namespace std::chrono;
        high_resolution_clock::time_point tend = high_resolution_clock::now();
        double dt = (duration_cast<duration<double>>(tend - tstart)).count();
        dt *= 1000.0;
        GlobalLog("%s: %f ms\n", label, dt);
    }

private:
    std::chrono::high_resolution_clock::time_point tstart;
};

template<int N> class MultiStopwatch {
public:
    MultiStopwatch()
    {
        for(int i = 0; i < N; ++i) {
            counts[i]       = 0;
            accumulators[i] = 0.0;
            maxes[i]        = 0.0;
            mins[i]         = 1e100;
        }
    }
    inline void tick(int i) { tstart[i] = std::chrono::high_resolution_clock::now(); }
    inline void tock(int i)
    {
        using namespace std::chrono;
        high_resolution_clock::time_point tend = high_resolution_clock::now();
        double dt = (duration_cast<duration<double>>(tend - tstart[i])).count();
        dt *= 1000.0;
        accumulators[i] += dt;
        if(maxes[i] < dt) maxes[i] = dt;
        if(mins[i] > dt) mins[i] = dt;
        ++counts[i];
    }
    inline void print(const char *label, const char *names[N])
    {
        std::stringstream ss;
        ss << label << ": " << std::setw(5);
        for(int i = 0; i < N; ++i) {
            ss << names[i] << " = (" << counts[i] << "/" << mins[i] << "/"
               << accumulators[i] << "/" << maxes[i] << ") ms, ";
        }
        GlobalLog("%s\n", ss.str().c_str());
    }

private:
    std::chrono::high_resolution_clock::time_point tstart[N];
    int counts[N];
    double accumulators[N];
    double maxes[N];
    double mins[N];
};

} // namespace bhc
