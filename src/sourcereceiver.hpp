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

template<bool O3D, bool R3D> inline void DefaultfreqVec(bhcParams<O3D, R3D> &params)
{
    params.freqinfo->Nfreq = 1;
    trackallocate(
        params, "default source frequencies", freqinfo->freqVec, params.freqinfo->Nfreq);
    freqinfo->freqVec[0] = freqinfo->freq0;
}

/**
 * Optionally reads a vector of source frequencies for a broadband run
 * If the broadband option is not selected, then the input freq (a scalar) is stored in
 * the frequency vector
 */
template<bool O3D, bool R3D> inline void ReadfreqVec(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile  = GetInternal(params)->PRTFile;
    FreqInfo *freqinfo     = params.freqinfo;
    params.freqinfo->Nfreq = 1;

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

template<bool O3D, bool R3D> inline void DefaultSzRz(bhcParams<O3D, R3D> &params)
{
    params.Pos->NSz = 1;
    trackallocate(params, "default source z-coordinates", params.Pos->Sz, 1);
    params.Pos->Sz[0] = RL(567.8);

    params.Pos->NRz = 11;
    trackallocate(
        params, "default receiver z-coordinates", params.Pos->Rz, params.Pos->NRz);
    for(int32_t i = 0; i < params.Pos->NRz; ++i) {
        params.Pos->Rz[i] = RL(500.0) * (real)i;
    }
}

/**
 * Reads source and receiver z-coordinates (depths)
 * zMin, zMax: limits for those depths;
 *     sources and receivers are shifted to be within those limits
 */
template<bool O3D, bool R3D> inline void ReadSzRz(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    Position *Pos         = params.Pos;
    real zMin             = params.Bdry->Top.hs.Depth;
    real zMax             = params.Bdry->Bot.hs.Depth;

    // bool monotonic; //LP: monotonic is a function, this is a name clash

    ReadVector(params, Pos->NSz, Pos->Sz, "Source   z-coordinates, Sz", "m", ENVFile);
    ReadVector(params, Pos->NRz, Pos->Rz, "Receiver z-coordinates, Rz", "m", ENVFile);

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

template<bool O3D, bool R3D> inline void DefaultRcvrRanges(bhcParams<O3D, R3D> &params)
{
    params.Pos->NRr = 10;
    trackallocate(
        params, "default receiver r-coordinates", params.Pos->Rr, params.Pos->NRr);
    for(int32_t i = 0; i < params.Pos->NRr; ++i) {
        params.Pos->Rr[i] = RL(5000.0) * (real)(i + 1);
    }
    Pos->Delta_r = Pos->Rr[Pos->NRr - 1] - Pos->Rr[Pos->NRr - 2];
}

template<bool O3D, bool R3D> inline void ReadRcvrRanges(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    Position *Pos = params.Pos;

    ReadVector(params, Pos->NRr, Pos->Rr, "Receiver r-coordinates, Rr", "km", ENVFile);

    // calculate range spacing TODO move to preprocessing
    Pos->Delta_r = FL(0.0);
    if(Pos->NRr != 1) Pos->Delta_r = Pos->Rr[Pos->NRr - 1] - Pos->Rr[Pos->NRr - 2];

    if(!monotonic(Pos->Rr, Pos->NRr)) {
        EXTERR("ReadRcvrRanges: Receiver ranges are not monotonically increasing");
    }
}

template<bool O3D, bool R3D> inline void DefaultRcvrBearings(bhcParams<O3D, R3D> &params)
{
    params.Pos->Ntheta = 5;
    trackallocate(
        params, "default receiver bearings", params.Pos->theta, params.Pos->Ntheta);
    trackallocate(params, "receiver bearing sin/cos table", Pos->t_rcvr, Pos->Ntheta);
    for(int32_t i = 0; i < params.Pos->Ntheta; ++i) {
        params.Pos->theta[i] = RL(72.0) * (real)i;
    }
    Pos->Delta_theta = Pos->theta[Pos->Ntheta - 1] - Pos->theta[Pos->Ntheta - 2];
}

template<bool O3D, bool R3D> inline void ReadRcvrBearings(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    Position *Pos = params.Pos;

    ReadVector(
        params, Pos->Ntheta, Pos->theta, "Receiver bearings, theta", "degrees", ENVFile);
    // TODO move this also to preprocessing
    trackallocate(params, "receiver bearing sin/cos table", Pos->t_rcvr, Pos->Ntheta);

    CheckFix360Sweep(Pos->theta, Pos->Ntheta);

    // calculate angular spacing TODO move to preprocessing
    Pos->Delta_theta = FL(0.0);
    if(Pos->Ntheta != 1)
        Pos->Delta_theta = Pos->theta[Pos->Ntheta - 1] - Pos->theta[Pos->Ntheta - 2];

    if(!monotonic(Pos->theta, Pos->Ntheta)) {
        EXTERR("ReadRcvrBearings: Receiver bearings are not monotonically increasing");
    }
}

} // namespace bhc
