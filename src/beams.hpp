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
#include "common.hpp"
#include "runtype.hpp"

namespace bhc {

template<bool O3D, bool R3D> inline void ReadPat(bhcParams<O3D, R3D> &params)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    BeamInfo *beaminfo    = params.beaminfo;

    params.beaminfo->SBPFlag = params.Beam->RunType[2];

    if(beaminfo->SBPFlag == '*') {
        PRTFile << "\n______________________________\nUsing source beam pattern file\n";

        LDIFile SBPFile(GetInternal(params), GetInternal(params)->FileRoot, ".sbp");
        if(!SBPFile.Good()) {
            PRTFile << "SBPFile = " << GetInternal(params)->FileRoot << ".sbp\n";
            EXTERR("BELLHOP-ReadPat: Unable to open source beampattern file");
        }

        LIST(SBPFile);
        SBPFile.Read(beaminfo->NSBPPts);
        PRTFile << "Number of source beam pattern points " << beaminfo->NSBPPts << "\n";

        trackallocate(
            params, "source beam pattern", beaminfo->SrcBmPat, beaminfo->NSBPPts * 2);

        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);

        for(int32_t i = 0; i < beaminfo->NSBPPts; ++i) {
            LIST(SBPFile);
            SBPFile.Read(&beaminfo->SrcBmPat[2 * i], 2);
            PRTFile << beaminfo->SrcBmPat[2 * i] << " " << beaminfo->SrcBmPat[2 * i + 1]
                    << "\n";
        }
    } else {
        beaminfo->NSBPPts = 2;
        trackallocate(
            params, "default/trivial source beam pattern", beaminfo->SrcBmPat, 2 * 2);
        beaminfo->SrcBmPat[0 * 2 + 0] = FL(-180.0);
        beaminfo->SrcBmPat[0 * 2 + 1] = FL(0.0);
        beaminfo->SrcBmPat[1 * 2 + 0] = FL(180.0);
        beaminfo->SrcBmPat[1 * 2 + 1] = FL(0.0);
    }

    if(!monotonic(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0)) {
        EXTERR("BELLHOP-ReadPat: Source beam pattern angles are not monotonic");
    }

    // convert dB to linear scale
    for(int32_t i = 0; i < beaminfo->NSBPPts; ++i)
        beaminfo->SrcBmPat[i * 2 + 1]
            = STD::pow(FL(10.0), beaminfo->SrcBmPat[i * 2 + 1] / FL(20.0));
}

template<bool O3D> inline HOST_DEVICE VEC23<O3D> BeamBoxCenter(const VEC23<O3D> &xs)
{
    VEC23<O3D> ret = xs;
    // box is centered at z=0
    DEP(ret) = RL(0.0);
    return ret;
}

template<bool O3D, int DIM> inline HOST_DEVICE bool IsOutsideBeamBoxDim(
    const VEC23<O3D> &x, const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs)
{
    static_assert(DIM >= 0 && DIM <= ZDIM<O3D>(), "Invalid use of IsOutsideBoxDim!");
    // LP: In 2D, source range is always 0.
    return STD::abs(x[DIM] - BeamBoxCenter<O3D>(xs)[DIM]) >= Beam->Box[DIM];
}

} // namespace bhc
