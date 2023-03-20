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
#include "paramsmodule.hpp"

namespace bhc {

/**
 * Read the RunType variable and echo with explanatory information to the print file
 */
template<bool O3D, bool R3D> class RunType {
public:
    RunType() {}
    virtual ~RunType() {}

    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        // RunType, infl/beam type, ignored, point source, rectilinear grid, dim, ignored
        memcpy(params.Beam->RunType, "CG RR  ", 7);
        if constexpr(R3D) {
            params.Beam->RunType[5] = '3';
        } else if constexpr(O3D) {
            params.Beam->RunType[5] = '2';
        }
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.Beam->RunType, 7);
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        switch(params.Beam->RunType[0]) {
        case 'R': break;
        case 'E': break;
        case 'I': break;
        case 'S': break;
        case 'C': break;
        case 'A': break;
        case 'a': break;
        default: EXTERR("ReadEnvironment: Unknown RunType selected");
        }

        switch(params.Beam->RunType[1]) {
        case 'C': break;
        case 'R': break;
        case 'S': break;
        case 'b': break;
        case 'B': break;
        case 'g': break;
        default: params.Beam->RunType[1] = 'G';
        }

        if(params.Beam->RunType[3] != 'X') params.Beam->RunType[3] = 'R';

        if(params.Beam->RunType[4] != 'I') params.Beam->RunType[4] = 'R';

        switch(params.Beam->RunType[5]) {
        case '3':
            if constexpr(!R3D) {
                EXTERR(
                    "Environment file specifies 3D, but you are running " BHC_PROGRAMNAME
                    " in 2D or Nx2D mode");
            }
            break;
        case '2':
            if constexpr(R3D) {
                EXTERR("Environment file specifies Nx2D or possibly 2D, but you are "
                       "running " BHC_PROGRAMNAME " in 3D mode");
            } else if constexpr(!O3D) {
                EXTWARN("Environment file specifies dimensionality 2, which usually "
                        "means Nx2D, but you are running " BHC_PROGRAMNAME " in 2D mode");
            }
            break;
        case ' ':
            if constexpr(R3D) {
                EXTERR("Environment file specifies 2D or possibly Nx2D, but you are "
                       "running " BHC_PROGRAMNAME " in 3D mode");
            } else if constexpr(O3D) {
                EXTWARN("Environment file specifies dimensionality ' ', which usually "
                        "means 2D, but you are running " BHC_PROGRAMNAME " in Nx2D mode");
            }
            break;
        default:
            EXTERR(
                "Unknown dimensionality %c in environment file", params.Beam->RunType[5]);
        }
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        PRTFile << "\n";

        switch(params.Beam->RunType[0]) {
        case 'R': PRTFile << "Ray trace run\n"; break;
        case 'E': PRTFile << "Eigenray trace run\n"; break;
        case 'I': PRTFile << "Incoherent TL calculation\n"; break;
        case 'S': PRTFile << "Semi-coherent TL calculation\n"; break;
        case 'C': PRTFile << "Coherent TL calculation\n"; break;
        case 'A': PRTFile << "Arrivals calculation, ASCII  file output\n"; break;
        case 'a': PRTFile << "Arrivals calculation, binary file output\n"; break;
        }

        switch(params.Beam->RunType[1]) {
        case 'C': PRTFile << "Cartesian beams\n"; break;
        case 'R': PRTFile << "Ray centered beams\n"; break;
        case 'S': PRTFile << "Simple gaussian beams\n"; break;
        case 'b':
            PRTFile << "Geometric gaussian beams in ray-centered coordinates\n";
            break;
        case 'B': PRTFile << "Geometric gaussian beams in Cartesian coordinates\n"; break;
        case 'g': PRTFile << "Geometric hat beams in ray-centered coordinates\n"; break;
        case 'G': PRTFile << "Geometric hat beams in Cartesian coordinates\n"; break;
        }

        switch(params.Beam->RunType[3]) {
        case 'X': PRTFile << "Line source (Cartesian coordinates)\n"; break;
        case 'R': PRTFile << "Point source (cylindrical coordinates)\n"; break;
        }

        switch(params.Beam->RunType[4]) {
        case 'I':
            PRTFile << "Irregular grid: Receivers at Rr[:] x Rz[:]\n";
            if(params.Pos->NRz != params.Pos->NRr)
                EXTWARN("ReadEnvironment: Irregular grid option selected with NRz not "
                        "equal to Nr");
            break;
        case 'R':
            PRTFile << "Rectilinear receiver grid: Receivers at Rr[:] x Rz[:]\n";
            break;
        }

        switch(params.Beam->RunType[5]) {
        case '3': PRTFile << "3D calculation\n"; break;
        case '2':
            PRTFile << "N x 2D calculation (neglects horizontal refraction)\n";
            break;
            // LP: No message for 2D mode.
        }
    }
};

} // namespace bhc