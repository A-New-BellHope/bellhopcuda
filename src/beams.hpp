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

inline void ReadPat(std::string FileRoot, PrintFileEmu &PRTFile, BeamInfo *beaminfo)
{
    if(beaminfo->SBPFlag == '*') {
        PRTFile << "\n______________________________\nUsing source beam pattern file\n";

        LDIFile SBPFile(GetInternal(params), FileRoot, ".sbp");
        if(!SBPFile.Good()) {
            PRTFile << "SBPFile = " << FileRoot << ".sbp\n";
            EXTERR("BELLHOP-ReadPat: Unable to open source beampattern file");
        }

        LIST(SBPFile);
        SBPFile.Read(beaminfo->NSBPPts);
        PRTFile << "Number of source beam pattern points " << beaminfo->NSBPPts << "\n";

        checkallocate(beaminfo->SrcBmPat, beaminfo->NSBPPts * 2);

        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);

        for(int32_t i = 0; i < beaminfo->NSBPPts; ++i) {
            LIST(SBPFile);
            SBPFile.Read(&beaminfo->SrcBmPat[2 * i], 2);
            PRTFile << beaminfo->SrcBmPat[2 * i] << " " << beaminfo->SrcBmPat[2 * i + 1]
                    << "\n";
        }
    } else {
        beaminfo->NSBPPts = 2;
        checkallocate(beaminfo->SrcBmPat, 2 * 2);
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
    return STD::abs(x[DIM] - BeamBoxCenter<O3D>(xs)[DIM]) > Beam->Box[DIM];
}

/**
 * Limits for tracing beams
 */
template<bool O3D> inline void ReadBeamInfo(
    LDIFile &ENVFile, PrintFileEmu &PRTFile, BeamStructure<O3D> *Beam,
    const BdryType *Bdry)
{
    if constexpr(O3D) {
        LIST(ENVFile);
        ENVFile.Read(Beam->deltas);
        ENVFile.Read(Beam->Box.x);
        ENVFile.Read(Beam->Box.y);
        ENVFile.Read(Beam->Box.z);
        Beam->Box.x *= FL(1000.0); // convert km to m
        Beam->Box.y *= FL(1000.0); // convert km to m

        if(Beam->deltas == FL(0.0))
            Beam->deltas = (Bdry->Bot.hs.Depth - Bdry->Top.hs.Depth)
                / FL(10.0); // Automatic step size selection
    } else {
        LIST(ENVFile);
        ENVFile.Read(Beam->deltas);
        ENVFile.Read(Beam->Box.y);
        ENVFile.Read(Beam->Box.x);
    }

    PRTFile << std::setprecision(4);
    PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam->deltas
            << " m\n\n";
    if constexpr(O3D) {
        PRTFile << "Maximum ray x-range, Box.x  = " << std::setw(11) << Beam->Box.x
                << " m\n";
        PRTFile << "Maximum ray y-range, Box.y  = " << std::setw(11) << Beam->Box.y
                << " m\n";
        PRTFile << "Maximum ray z-range, Box.z  = " << std::setw(11) << Beam->Box.z
                << " m\n";
    } else {
        PRTFile << "Maximum ray depth, Box.y  = " << std::setw(11) << Beam->Box.y
                << " m\n";
        PRTFile << "Maximum ray range, Box.x  = " << std::setw(11) << Beam->Box.x
                << "km\n";

        Beam->Box.x *= FL(1000.0); // convert km to m
    }

    // *** Beam characteristics ***

    Beam->Type[3] = Beam->RunType[6]; // selects beam shift option

    if(Beam->Type[3] == 'S') {
        PRTFile << "Beam shift in effect\n";
    } else {
        PRTFile << "No beam shift in effect\n";
    }

    if(!IsRayRun(Beam)) { // no worry about the beam type if this is a ray trace run

        // Curvature change can cause overflow in grazing case
        // Suppress by setting BeamType( 3 : 3 ) = 'Z'

        Beam->Type[0] = Beam->RunType[1];
        if(IsGeometricInfl(Beam) || IsSGBInfl(Beam)) {
            NULLSTATEMENT;
        } else if(IsCervenyInfl(Beam)) {
            LIST(ENVFile);
            ENVFile.Read(&Beam->Type[1], 2);
            ENVFile.Read(Beam->epsMultiplier);
            ENVFile.Read(Beam->rLoop);
            PRTFile << "\n\nType of beam = " << Beam->Type[0] << "\n";
            switch(Beam->Type[2]) {
            case 'D': PRTFile << "Curvature doubling invoked\n"; break;
            case 'Z': PRTFile << "Curvature zeroing invoked\n"; break;
            case 'S': PRTFile << "Standard curvature condition\n"; break;
            default:
                EXTERR(
                    "ReadEnvironment: Unknown curvature condition '%c'", Beam->Type[2]);
            }

            PRTFile << "Epsilon multiplier " << Beam->epsMultiplier << "\n";
            PRTFile << "Range for choosing beam width " << Beam->rLoop << "\n";

            // Images, windows
            // LP: These values are not initialized if not written in the file,
            // and Component is not always written in the test env files.
            Beam->Nimage      = 1;
            Beam->iBeamWindow = 4;
            Beam->Component   = 'P';
            LIST(ENVFile);
            ENVFile.Read(Beam->Nimage);
            ENVFile.Read(Beam->iBeamWindow);
            ENVFile.Read(Beam->Component);
            PRTFile << "\nNumber of images, Nimage  = " << Beam->Nimage << "\n";
            PRTFile << "Beam windowing parameter  = " << Beam->iBeamWindow << "\n";
            PRTFile << "Component                 = " << Beam->Component << "\n";
        } else {
            EXTERR(
                "ReadEnvironment: Unknown beam type (second letter of run type == '%c')",
                Beam->RunType[1]);
        }

        PRTFile << "\n";
    }
}

} // namespace bhc
