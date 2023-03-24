/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
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
#include "../common_setup.hpp"
#include "field.hpp"

namespace bhc { namespace mode {

template<bool O3D, bool R3D> void PostProcessArrivals(
    const bhcParams<O3D, R3D> &params, ArrInfo *arrinfo);
extern template void PostProcessArrivals<false, false>(
    const bhcParams<false, false> &params, ArrInfo *arrinfo);
extern template void PostProcessArrivals<true, false>(
    const bhcParams<true, false> &params, ArrInfo *arrinfo);
extern template void PostProcessArrivals<true, true>(
    const bhcParams<true, true> &params, ArrInfo *arrinfo);

template<bool O3D, bool R3D> void WriteOutArrivals(
    const bhcParams<O3D, R3D> &params, const ArrInfo *arrinfo);
extern template void WriteOutArrivals<false, false>(
    const bhcParams<false, false> &params, const ArrInfo *arrinfo);
extern template void WriteOutArrivals<true, false>(
    const bhcParams<true, false> &params, const ArrInfo *arrinfo);
extern template void WriteOutArrivals<true, true>(
    const bhcParams<true, true> &params, const ArrInfo *arrinfo);

template<bool O3D, bool R3D> class Arr : public Field {
public:
    Arr() {}
    virtual ~Arr() {}

    virtual void Init(bhcOutputs<O3D, R3D> &outputs) const
    {
        outputs.arrinfo->Arr           = nullptr;
        outputs.arrinfo->NArr          = nullptr;
        outputs.arrinfo->MaxNPerSource = nullptr;
        outputs.arrinfo->MaxNArr       = 1;
    }

    virtual void Preprocess(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        Field::Preprocess(params);

        trackdeallocate(params, arrinfo->Arr);
        trackdeallocate(params, arrinfo->NArr);
        trackdeallocate(params, arrinfo->MaxNPerSource);
        arrinfo->AllowMerging = GetInternal(params)->numThreads == 1;
        size_t nSrcs          = params.Pos->NSx * params.Pos->NSy * params.Pos->NSz;
        size_t nSrcsRcvrs     = nSrcs * params.Pos->Ntheta * params.Pos->NRr
            * params.Pos->NRz_per_range;
        int64_t remainingMemory = GetInternal(params)->maxMemory
            - GetInternal(params)->usedMemory;
        remainingMemory -= nSrcsRcvrs * sizeof(int32_t);
        remainingMemory -= nSrcs * sizeof(int32_t);
        remainingMemory -= 32 * 3; // Possible padding used for the three arrays
        remainingMemory  = std::max(remainingMemory, (int64_t)0);
        arrinfo->MaxNArr = remainingMemory / (nSrcsRcvrs * sizeof(Arrival));
        if(arrinfo->MaxNArr == 0) {
            EXTERR("Insufficient memory to allocate arrivals");
        } else if(arrinfo->MaxNArr < 10) {
            EXTWARN(
                "Only enough memory to allocate up to %d arrivals per receiver",
                arrinfo->MaxNArr);
        }
        GetInternal(params)->PRTFile << "\n( Maximum # of arrivals = " << arrinfo->MaxNArr
                                     << " )\n";
        trackallocate(
            params, "arrivals", arrinfo->Arr, nSrcsRcvrs * (size_t)arrinfo->MaxNArr);
        trackallocate(params, "arrivals", arrinfo->NArr, nSrcsRcvrs);
        trackallocate(params, "arrivals", arrinfo->MaxNPerSource, nSrcs);
        memset(arrinfo->Arr, 0, nSrcsRcvrs * (size_t)arrinfo->MaxNArr * sizeof(Arrival));
        memset(arrinfo->NArr, 0, nSrcsRcvrs * sizeof(int32_t));
        // MaxNPerSource does not have to be initialized
    }

    virtual void Postprocess(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        PostProcessArrivals<O3D, R3D>(params, outputs.arrinfo);
    }

    virtual void Writeout(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        WriteOutArrivals<O3D, R3D>(params, outputs.arrinfo);
    }

    virtual void Finalize(bhcOutputs<O3D, R3D> &outputs) const
    {
        trackdeallocate(params, outputs.arrinfo->Arr);
        trackdeallocate(params, outputs.arrinfo->NArr);
        trackdeallocate(params, outputs.arrinfo->MaxNPerSource);
    }
};

}} // namespace bhc::mode
