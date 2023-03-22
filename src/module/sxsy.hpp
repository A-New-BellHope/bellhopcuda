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

namespace bhc { namespace module {

/**
 * source x-y coordinates
 */
template<bool O3D, bool R3D> class SxSy : public ParamsModule {
public:
    SxSy() {}
    virtual ~SxSy() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const override
    {
        params.Pos->Sx = nullptr;
        params.Pos->Sy = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const override
    {
        params.Pos->SxSyInKm = true;
        params.Pos->NSx      = 1;
        params.Pos->NSy      = 1;
    }
    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        trackallocate(params, "default/trivial source x-coordinates", params.Pos->Sx, 1);
        trackallocate(params, "default/trivial source y-coordinates", params.Pos->Sy, 1);
        params.Pos->Sx[0] = FL(0.0); // dummy x-coordinate
        params.Pos->Sy[0] = FL(0.0); // dummy y-coordinate
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        if constexpr(O3D) {
            ReadVector2(params, params.Pos->NSx, params.Pos->Sx, ENVFile);
            ReadVector2(params, params.Pos->NSy, params.Pos->Sy, ENVFile);
        } else {
            Default(params);
        }
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        ValidateVector2(
            params, params.Pos->NSx, params.Pos->Sx, "Source   x-coordinates, Sx");
        ValidateVector2(
            params, params.Pos->NSy, params.Pos->Sy, "Source   y-coordinates, Sy");
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        Preprocess(params);
        EchoVector2(
            params, params.Pos->NSx, params.Pos->Sx, RL(0.001),
            "Source   x-coordinates, Sx", "km");
        EchoVector2(
            params, params.Pos->NSy, params.Pos->Sy, RL(0.001),
            "Source   y-coordinates, Sy", "km");
    }
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const override
    {
        if(!params.Pos->SxSyInKm) return;
        ToMeters2(params.Pos->NSx, params.Pos->Sx);
        ToMeters2(params.Pos->NSy, params.Pos->Sy);
        params.Pos->SxSyInKm = false;
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const override
    {
        trackdeallocate(params, params.Pos->Sx);
        trackdeallocate(params, params.Pos->Sy);
    }
};

}} // namespace bhc::module
