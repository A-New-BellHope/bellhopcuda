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
#include "run.hpp"

#include <thread>
#include <vector>

namespace bhc {

template<bool O3D, bool R3D> void FieldModesWorker(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    try {
        while(true) {
            int32_t job = sharedJobID++;
            RayInitInfo rinit;
            if(!GetJobIndices<O3D>(rinit, job, params.Pos, params.Angles)) break;

            MainFieldModes<O3D, R3D>(
                rinit, outputs.uAllSources, params.Bdry, params.bdinfo, params.refl,
                params.ssp, params.Pos, params.Angles, params.freqinfo, params.Beam,
                params.beaminfo, outputs.eigen, outputs.arrinfo);
        }
    } catch(const std::exception &e) {
        std::lock_guard<std::mutex> lock(exceptionMutex);
        exceptionStr += std::string(e.what()) + "\n";
    }
}

template<bool O3D, bool R3D> void RunFieldModesImpl(
    bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs, uint32_t cores)
{
    std::vector<std::thread> threads;
    for(uint32_t i = 0; i < cores; ++i)
        threads.push_back(
            std::thread(FieldModesWorker<O3D, R3D>, std::ref(params), std::ref(outputs)));
    for(uint32_t i = 0; i < cores; ++i) threads[i].join();
}

#if BHC_ENABLE_2D
template void RunFieldModesImpl<false, false>(
    bhcParams<false, false> &params, bhcOutputs<false, false> &outputs, uint32_t cores);
#endif
#if BHC_ENABLE_NX2D
template void RunFieldModesImpl<true, false>(
    bhcParams<true, false> &params, bhcOutputs<true, false> &outputs, uint32_t cores);
#endif
#if BHC_ENABLE_3D
template void RunFieldModesImpl<true, true>(
    bhcParams<true, true> &params, bhcOutputs<true, true> &outputs, uint32_t cores);
#endif

} // namespace bhc
