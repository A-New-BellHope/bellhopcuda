/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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
template<bool O3D> class SzRz : public ParamsModule<O3D> {
public:
    SzRz() {}
    virtual ~SzRz() {}

    virtual void Init(bhcParams<O3D> &params) const override
    {
        params.Pos->Sz = nullptr;
        params.Pos->Rz = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D> &params) const override
    {
        params.Pos->NSz = 1;
        params.Pos->NRz = 1;
    }
    virtual void Default(bhcParams<O3D> &params) const override
    {
        trackallocate(params, "default source z-coordinates", params.Pos->Sz, 1);
        params.Pos->Sz[0] = FL(567.8);

        params.Pos->NRz = 11;
        trackallocate(
            params, "default receiver z-coordinates", params.Pos->Rz, params.Pos->NRz);
        for(int32_t i = 0; i < params.Pos->NRz; ++i) {
            params.Pos->Rz[i] = RL(500.0) * (real)i;
        }
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        ReadVector(params, params.Pos->Sz, params.Pos->NSz, ENVFile, DescriptionS);
        ReadVector(params, params.Pos->Rz, params.Pos->NRz, ENVFile, DescriptionR);
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        UnSubTab(ENVFile, params.Pos->Sz, params.Pos->NSz, "NSD", "SD(1:NSD) (m)");
        UnSubTab(ENVFile, params.Pos->Rz, params.Pos->NRz, "NRD", "RD(1:NRD) (m)");
    }
    void ExtSetupSz(bhcParams<O3D> &params, int32_t NSz) const
    {
        params.Pos->NSz = NSz;
        trackallocate(params, DescriptionS, params.Pos->Sz, params.Pos->NSz);
    }
    void ExtSetupRz(bhcParams<O3D> &params, int32_t NRz) const
    {
        params.Pos->NRz = NRz;
        trackallocate(params, DescriptionR, params.Pos->Rz, params.Pos->NRz);
    }
    virtual void Validate(bhcParams<O3D> &params) const override
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

        ValidateVector(params, Pos->Sz, Pos->NSz, DescriptionS, true);
        ValidateVector(params, Pos->Rz, Pos->NRz, DescriptionR, true);
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        EchoVectorWDescr(
            params, params.Pos->Sz, params.Pos->NSz, FL(1.0), DescriptionS, Units);
        EchoVectorWDescr(
            params, params.Pos->Rz, params.Pos->NRz, FL(1.0), DescriptionR, Units);
    }
    virtual void Preprocess(bhcParams<O3D> &params) const override
    {
        // irregular or rectilinear grid
        params.Pos->NRz_per_range = IsIrregularGrid(params.Beam) ? 1 : params.Pos->NRz;
    }
    virtual void Finalize(bhcParams<O3D> &params) const override
    {
        trackdeallocate(params, params.Pos->Sz);
        trackdeallocate(params, params.Pos->Rz);
    }

private:
    constexpr static const char *DescriptionS = "Source   z-coordinates, Sz";
    constexpr static const char *DescriptionR = "Receiver z-coordinates, Rz";
    constexpr static const char *Units        = "m";
};

}} // namespace bhc::module
