/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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

namespace bhc {

// #define BHC_USE_HIGH_PRIORITY_THREADS 1

void SetupThread();

inline int32_t ModifyNumThreads(int32_t numThreads)
{
    if(numThreads >= 1) return numThreads;
    numThreads = std::thread::hardware_concurrency();
    if(numThreads < 1) numThreads = 1;
#ifdef BHC_USE_HIGH_PRIORITY_THREADS
    if(numThreads > 1) --numThreads; // Leave one core unused to avoid locking up system
#endif
    return numThreads;
}

class Stopwatch {
public:
    Stopwatch(bhcInternal *internal_) : internal(internal_) {}
    inline void tick() { tstart = std::chrono::high_resolution_clock::now(); }
    inline void tock(const char *label)
    {
        using namespace std::chrono;
        high_resolution_clock::time_point tend = high_resolution_clock::now();
        double dt = (duration_cast<duration<double>>(tend - tstart)).count();
        dt *= 1000.0;
        ExternalWarning(internal, "%s: %f ms", label, dt);
    }

private:
    bhcInternal *internal;
    std::chrono::high_resolution_clock::time_point tstart;
};

template<int N> class MultiStopwatch {
public:
    MultiStopwatch(bhcInternal *internal_) : internal(internal_)
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
        ExternalWarning(internal, "%s", ss.str().c_str());
    }

private:
    bhcInternal *internal;
    std::chrono::high_resolution_clock::time_point tstart[N];
    int counts[N];
    double accumulators[N];
    double maxes[N];
    double mins[N];
};

} // namespace bhc
