/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP / BELLHOP3D underwater acoustics
simulator Copyright (C) 2021-2023 The Regents of the University of California Marine
Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu Based on BELLHOP /
BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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

template<bool O3D, bool R3D> void PostProcessTL(
    const bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs);
extern template void PostProcessTL<false, false>(
    const bhcParams<false> &params, bhcOutputs<false, false> &outputs);
extern template void PostProcessTL<true, false>(
    const bhcParams<true> &params, bhcOutputs<true, false> &outputs);
extern template void PostProcessTL<true, true>(
    const bhcParams<true> &params, bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> void WriteOutTL(
    const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs);
extern template void WriteOutTL<false, false>(
    const bhcParams<false> &params, const bhcOutputs<false, false> &outputs);
extern template void WriteOutTL<true, false>(
    const bhcParams<true> &params, const bhcOutputs<true, false> &outputs);
extern template void WriteOutTL<true, true>(
    const bhcParams<true> &params, const bhcOutputs<true, true> &outputs);

template<bool O3D, bool R3D> class TL : public Field<O3D, R3D> {
public:
    TL() {}
    virtual ~TL() {}

    virtual void Init(bhcOutputs<O3D, R3D> &outputs) const
    {
        outputs.uAllSources = nullptr;
    }

    virtual void Preprocess(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        Field<O3D, R3D>::Preprocess(params, outputs);

        trackdeallocate(params, outputs.uAllSources); // Free if previously run
        // for a TL calculation, allocate space for the pressure matrix
        const Position *Pos = params.Pos;
        size_t n            = (size_t)Pos->NSz * (size_t)Pos->NSx * (size_t)Pos->NSy
            * (size_t)Pos->Ntheta * (size_t)Pos->NRz_per_range * (size_t)Pos->NRr;
        trackallocate(params, "sound field / transmission loss", outputs.uAllSources, n);
        memset(outputs.uAllSources, 0, n * sizeof(cpxf));
    }

    virtual void Postprocess(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        PostProcessTL<O3D, R3D>(params, outputs);
    }

    virtual void Writeout(
        const bhcParams<O3D> &params, const bhcOutputs<O3D, R3D> &outputs) const
    {
        WriteOutTL<O3D, R3D>(params, outputs);
    }

    virtual void Finalize(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        trackdeallocate(params, outputs.uAllSources);
    }
};

}} // namespace bhc::mode
