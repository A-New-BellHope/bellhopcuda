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
#include "paramsmodule.hpp"

namespace bhc { namespace module {

/**
 * source and receiver z-coordinates (depths)
 */
template<bool O3D, bool R3D> class SzRz : public ParamsModule<O3D, R3D> {
public:
    SzRz() {}
    virtual ~SzRz() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        params.Pos->Sz = nullptr;
        params.Pos->Rz = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const
    {
        params.Pos->NSz = 1;
        params.Pos->NRz = 1;
    }
    virtual void Default(bhcParams<O3D, R3D> &params) const
    {
        trackallocate(params, "default source z-coordinates", params.Pos->Sz, 1);
        params.Pos->Sz[0] = RL(567.8);

        params.Pos->NRz = 11;
        trackallocate(
            params, "default receiver z-coordinates", params.Pos->Rz, params.Pos->NRz);
        for(int32_t i = 0; i < params.Pos->NRz; ++i) {
            params.Pos->Rz[i] = RL(500.0) * (real)i;
        }
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        ReadVector(params, params.Pos->NSz, params.Pos->Sz, ENVFile, DescriptionS);
        ReadVector(params, params.Pos->NRz, params.Pos->Rz, ENVFile, DescriptionR);
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        Position *Pos         = params.Pos;

        // zMin, zMax: limits for those depths;
        //     sources and receivers are shifted to be within those limits
        real zMin = params.Bdry->Top.hs.Depth;
        real zMax = params.Bdry->Bot.hs.Depth;

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
            PRTFile
                << "Warning in ReadSzRz : Source above or too near the top bdry has been "
                   "moved down\n";
        if(botbdry)
            PRTFile
                << "Warning in ReadSzRz : Source below or too near the bottom bdry has "
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
            PRTFile
                << "Warning in ReadSzRz : Receiver above or too near the top bdry has "
                   "been moved down\n";
        if(botbdry)
            PRTFile
                << "Warning in ReadSzRz : Receiver below or too near the bottom bdry has "
                   "been moved up\n";

        ValidateVector(params, Pos->NSz, Pos->Sz, DescriptionS);
        ValidateVector(params, Pos->NRz, Pos->Rz, DescriptionR);
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        EchoVectorWDescr(params, params.Pos->NSz, params.Pos->Sz, DescriptionS, Units);
        EchoVectorWDescr(params, params.Pos->NRz, params.Pos->Rz, DescriptionR, Units);
    }
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const
    {
        // irregular or rectilinear grid
        params.Pos->NRz_per_range = IsIrregularGrid(params.Beam) ? 1 : params.Pos->NRz;
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        trackdeallocate(params, params.Pos->Sz);
        trackdeallocate(params, params.Pos->Rz);
    }

private:
    constexpr const char *DescriptionS = "Source   z-coordinates, Sz";
    constexpr const char *DescriptionR = "Receiver z-coordinates, Rz";
    constexpr const char *Units        = "m";
};

}} // namespace bhc::module
