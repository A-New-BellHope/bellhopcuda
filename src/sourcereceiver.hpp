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
 * Optionally reads a vector of source frequencies for a broadband run
 * If the broadband option is not selected, then the input freq (a scalar) is stored in
 * the frequency vector
 */
template<bool O3D, bool R3D> inline void ReadfreqVec(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    FreqInfo *freqinfo    = params.freqinfo;

    if(params.Bdry->Top.hs.Opt[5] == 'B') {
        LIST(ENVFile);
        ENVFile.Read(freqinfo->Nfreq);
        PRTFile << "_____________________________________________________________________"
                   "_____\n\n\n";
        PRTFile << "   Number of frequencies = " << freqinfo->Nfreq << "\n";
        if(freqinfo->Nfreq <= 0) { EXTERR("Number of frequencies must be positive"); }
    }

    trackallocate(
        params, "source frequencies", freqinfo->freqVec, bhc::max(3, freqinfo->Nfreq));

    if(params.Bdry->Top.hs.Opt[5] == 'B') {
        PRTFile << "   Frequencies (Hz)\n";
        freqinfo->freqVec[1] = FL(-999.9);
        freqinfo->freqVec[2] = FL(-999.9);
        LIST(ENVFile);
        ENVFile.Read(freqinfo->freqVec, freqinfo->Nfreq);
        SubTab(freqinfo->freqVec, freqinfo->Nfreq);
        EchoVector(freqinfo->freqVec, freqinfo->Nfreq, PRTFile);
    } else {
        freqinfo->freqVec[0] = freqinfo->freq0;
    }
}

/**
 * Read a vector x
 * Description is something like 'receiver ranges'
 * Units       is something like 'km'
 */
template<bool O3D, bool R3D, typename REAL> inline void ReadVector(
    bhcParams<O3D, R3D> &params, int32_t &Nx, REAL *&x, std::string Description,
    std::string Units, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    PRTFile << "\n_______________________________________________________________________"
               "___\n\n";
    LIST(ENVFile);
    ENVFile.Read(Nx);
    PRTFile << "   Number of " << Description << " = " << Nx << "\n";

    if(Nx <= 0) {
        EXTERR("ReadVector: Number of %s must be positive", Description.c_str());
    }

    trackallocate(params, Description.c_str(), x, bhc::max(3, Nx));

    PRTFile << "   " << Description << " (" << Units << ")\n";
    x[2] = FL(-999.9);
    LIST(ENVFile);
    ENVFile.Read(x, Nx);

    SubTab(x, Nx);
    Sort(x, Nx);
    EchoVector(x, Nx, PRTFile, 10, "   ");

    PRTFile << "\n";

    // Vectors in km should be converted to m for internal use
    trim(Units);
    if(Units.length() >= 2) {
        if(Units.substr(0, 2) == "km")
            for(int32_t i = 0; i < Nx; ++i) x[i] *= FL(1000.0);
    }
}

/**
 * Read source x-y coordinates
 *
 * ThreeD: flag indicating whether this is a 3D run
 */
template<bool O3D, bool R3D> inline void ReadSxSy(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    if constexpr(O3D) {
        ReadVector(
            params, params.Pos->NSx, params.Pos->Sx, "Source   x-coordinates, Sx", "km",
            ENVFile);
        ReadVector(
            params, params.Pos->NSy, params.Pos->Sy, "Source   y-coordinates, Sy", "km",
            ENVFile);
    } else {
        trackallocate(params, "default/trivial source x-coordinates", params.Pos->Sx, 1);
        trackallocate(params, "default/trivial source y-coordinates", params.Pos->Sy, 1);
        params.Pos->Sx[0] = FL(0.0);
        params.Pos->Sy[0] = FL(0.0);
    }
}

/**
 * Reads source and receiver z-coordinates (depths)
 * zMin, zMax: limits for those depths;
 *     sources and receivers are shifted to be within those limits
 */
template<bool O3D, bool R3D> inline void ReadSzRz(
    bhcParams<O3D, R3D> &params, real zMin, real zMax, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    Position *Pos         = params.Pos;

    // bool monotonic; //LP: monotonic is a function, this is a name clash

    ReadVector(params, Pos->NSz, Pos->Sz, "Source   z-coordinates, Sz", "m", ENVFile);
    ReadVector(params, Pos->NRz, Pos->Rz, "Receiver z-coordinates, Rz", "m", ENVFile);

    trackallocate(params, "source depth auxiliary data", Pos->ws, Pos->NSz);
    trackallocate(params, "source depth auxiliary data", Pos->iSz, Pos->NSz);
    trackallocate(params, "receiver depth auxiliary data", Pos->wr, Pos->NRz);
    trackallocate(params, "receiver depth auxiliary data", Pos->iRz, Pos->NRz);

    // *** Check for Sz/Rz in upper or lower halfspace ***

    bool topbdry = false, botbdry = false;
    for(int32_t i = 0; i < Pos->NSz; ++i) {
        if(Pos->Sz[i] < zMin) {
            topbdry    = true;
            Pos->Sz[i] = zMin;
        }
        if(Pos->Sz[i] > zMax) {
            botbdry    = true;
            Pos->Sz[i] = zMax;
        }
    }
    if(topbdry)
        PRTFile << "Warning in ReadSzRz : Source above or too near the top bdry has been "
                   "moved down\n";
    if(botbdry)
        PRTFile << "Warning in ReadSzRz : Source below or too near the bottom bdry has "
                   "been moved up\n";

    topbdry = false;
    botbdry = false;
    for(int32_t i = 0; i < Pos->NRz; ++i) {
        if(Pos->Rz[i] < zMin) {
            topbdry    = true;
            Pos->Rz[i] = zMin;
        }
        if(Pos->Rz[i] > zMax) {
            botbdry    = true;
            Pos->Rz[i] = zMax;
        }
    }
    if(topbdry)
        PRTFile << "Warning in ReadSzRz : Receiver above or too near the top bdry has "
                   "been moved down\n";
    if(botbdry)
        PRTFile << "Warning in ReadSzRz : Receiver below or too near the bottom bdry has "
                   "been moved up\n";

    /*
    if(!monotonic(Pos->sz, Pos->NSz)){
        EXTERR("SzRzRMod: Source depths are not monotonically increasing");
    }
    if(!monotonic(Pos->rz, Pos->NRz)){
        EXTERR("SzRzRMod: Receiver depths are not monotonically increasing");
    }
    */
}

template<bool O3D, bool R3D> inline void ReadRcvrRanges(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    Position *Pos = params.Pos;

    ReadVector(params, Pos->NRr, Pos->Rr, "Receiver r-coordinates, Rr", "km", ENVFile);

    // calculate range spacing
    Pos->Delta_r = FL(0.0);
    if(Pos->NRr != 1) Pos->Delta_r = Pos->Rr[Pos->NRr - 1] - Pos->Rr[Pos->NRr - 2];

    if(!monotonic(Pos->Rr, Pos->NRr)) {
        EXTERR("ReadRcvrRanges: Receiver ranges are not monotonically increasing");
    }
}

template<bool O3D, bool R3D> inline void ReadRcvrBearings(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    Position *Pos = params.Pos;

    ReadVector(
        params, Pos->Ntheta, Pos->theta, "Receiver bearings, theta", "degrees", ENVFile);
    trackallocate(params, "receiver bearing sin/cos table", Pos->t_rcvr, Pos->Ntheta);

    CheckFix360Sweep(Pos->theta, Pos->Ntheta);

    // calculate angular spacing
    Pos->Delta_theta = FL(0.0);
    if(Pos->Ntheta != 1)
        Pos->Delta_theta = Pos->theta[Pos->Ntheta - 1] - Pos->theta[Pos->Ntheta - 2];

    if(!monotonic(Pos->theta, Pos->Ntheta)) {
        EXTERR("ReadRcvrBearings: Receiver bearings are not monotonically increasing");
    }
}

} // namespace bhc
