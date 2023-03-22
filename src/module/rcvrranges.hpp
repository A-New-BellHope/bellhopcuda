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

template<bool O3D, bool R3D> class RcvrRanges : public ParamsModule {
public:
    RcvrRanges() {}
    virtual ~RcvrRanges() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const override
    {
        params.Pos->Rr = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const override
    {
        params.Pos->NRr = 1;
    }
    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        params.Pos->NRr = 10;
        trackallocate(
            params, "default receiver r-coordinates", params.Pos->Rr, params.Pos->NRr);
        for(int32_t i = 0; i < params.Pos->NRr; ++i) {
            params.Pos->Rr[i] = RL(5000.0) * (real)(i + 1);
        }
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        ReadVector2(params, params.Pos->NRr, params.Pos->Rr, ENVFile);
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        ValidateVector2(params, params.Pos->NRr, params.Pos->Rr, "Receiver ranges");
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        Preprocess(params);
        EchoVector2(
            params, params.Pos->NRr, params.Pos->Rr, RL(0.001),
            "Receiver r-coordinates, Rr", "km");
    }
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const override
    {
        if(params.Pos->RrInKm) {
            ToMeters2(params.Pos->NRr, params.Pos->Rr);
            params.Pos->RrInKm = false;
        }

        // calculate range spacing
        Pos->Delta_r = FL(0.0);
        if(Pos->NRr >= 2) Pos->Delta_r = Pos->Rr[Pos->NRr - 1] - Pos->Rr[Pos->NRr - 2];
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const override
    {
        trackdeallocate(params, params.Pos->Rr);
    }
};

}} // namespace bhc::module
