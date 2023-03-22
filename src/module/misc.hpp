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

namespace bhc {

/**
 * Read a vector x
 */
template<bool O3D, bool R3D, typename REAL> inline void ReadVector2(
    bhcParams<O3D, R3D> &params, int32_t &Nx, REAL *&x, LDIFile &ENVFile)
{
    LIST(ENVFile);
    ENVFile.Read(Nx);
    trackallocate(params, Description.c_str(), x, bhc::max(3, Nx));
    x[1] = FL(-999.9);
    x[2] = FL(-999.9);
    LIST(ENVFile);
    ENVFile.Read(x, Nx);
    SubTab(x, Nx);
    Sort(x, Nx);
}

template<bool O3D, bool R3D, typename REAL> inline void ValidateVector2(
    bhcParams<O3D, R3D> &params, int32_t &Nx, REAL *&x, std::string Description)
{
    if(Nx <= 0) {
        EXTERR("ValidateVector2: Number of %s must be positive", Description.c_str());
    }
    if(!monotonic(x, Nx)) {
        EXTERR(
            "ValidateVector2: %s are not monotonically increasing", Description.c_str());
    }
}

/**
 * Description is something like 'receiver ranges'
 * Units       is something like 'km'
 */
template<bool O3D, bool R3D, typename REAL> inline void EchoVector2(
    bhcParams<O3D, R3D> &params, int32_t &Nx, REAL *&x, REAL multiplier,
    std::string Description, std::string Units)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    PRTFile << "\n_______________________________________________________________________"
               "___\n\n";
    PRTFile << "   Number of " << Description << " = " << Nx << "\n";
    PRTFile << "   " << Description << " (" << Units << ")\n";
    EchoVector(x, Nx, PRTFile, 10, "   ", multiplier);
    PRTFile << "\n";
}

template<typename REAL> inline void ToMeters2(int32_t &Nx, REAL *&x)
{
    for(int32_t i = 0; i < Nx; ++i) x[i] *= FL(1000.0);
}

} // namespace bhc
