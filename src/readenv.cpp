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
#include "readenv.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

namespace bhc {

template<bool O3D, bool R3D> void ReadEnvironment(bhcParams<O3D, R3D> &params)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    // Values only initialized once--reused from top to ssp, and ssp to bot
    HSInfo RecycledHS;
    RecycledHS.alphaR = FL(1500.0);
    RecycledHS.betaR  = FL(0.0);
    RecycledHS.alphaI = FL(0.0);
    RecycledHS.betaI  = FL(0.0);
    RecycledHS.rho    = FL(1.0);

    PRTFile << BHC_PROGRAMNAME << (R3D ? "3D" : O3D ? "Nx2D" : "") << "\n\n";

    // Open the environmental file
    LDIFile ENVFile(GetInternal(params), GetInternal(params)->FileRoot + ".env");
    if(!ENVFile.Good()) {
        PRTFile << "ENVFile = " << GetInternal(params)->FileRoot << ".env\n";
        EXTERR(BHC_PROGRAMNAME
               " - ReadEnvironment: Unable to open the environmental file");
    }

    ReadTitle<O3D, R3D>(params, ENVFile);
    ReadFreq0<O3D, R3D>(params, ENVFile);
    ReadNMedia<O3D, R3D>(params, ENVFile);
    ReadTopOpt<O3D, R3D>(params, ENVFile);
    BoundaryCondTop<O3D, R3D>(params, ENVFile, RecycledHS);
    ReadSSP<O3D, R3D>(params, ENVFile, RecycledHS);
    ReadBotOpt<O3D, R3D>(params, ENVFile);
    BoundaryCondBot<O3D, R3D>(params, ENVFile, RecycledHS);
    ReadSxSy<O3D, R3D>(params, ENVFile);
    ReadSzRz<O3D, R3D>(params, ENVFile);
    ReadRcvrRanges<O3D, R3D>(params, ENVFile);
    ReadRcvrBearings<O3D, R3D>(params, ENVFile);
    ReadfreqVec<O3D, R3D>(params, ENVFile);
    ReadRunType<O3D, R3D>(params, ENVFile);
    ReadRayAnglesElevation<O3D, R3D>(params, ENVFile);
    ReadRayAnglesBearing<O3D, R3D>(params, ENVFile);
    ReadBeamInfo<O3D, R3D>(params, ENVFile);
    ReadAltimetry<O3D, R3D>(params, ENVFile);
    ReadBathymetry<O3D, R3D>(params, ENVFile);
    ReadReflectionCoefficient<O3D, R3D>(params); //
    ReadPat<O3D, R3D>(params);                   //

    GetInternal(params)->PRTFile << "\n";
}

#if BHC_ENABLE_2D
template void ReadEnvironment<false, false>(
    bhcParams<false, false> &params, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_NX2D
template void ReadEnvironment<true, false>(
    bhcParams<true, false> &params, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_3D
template void ReadEnvironment<true, true>(
    bhcParams<true, true> &params, HSInfo &RecycledHS);
#endif

} // namespace bhc
