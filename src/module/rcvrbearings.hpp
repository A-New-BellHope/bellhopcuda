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
 * source x-y coordinates
 */
template<bool O3D, bool R3D> class RcvrBearings : public ParamsModule {
public:
    RcvrBearings() {}
    virtual ~RcvrBearings() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        params.Pos->theta  = nullptr;
        params.Pos->t_rcvr = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const { params.Pos->Ntheta = 1; }
    virtual void Default(bhcParams<O3D, R3D> &params) const
    {
        if constexpr(O3D) { params.Pos->Ntheta = 5; }
        trackallocate(
            params, "default receiver bearings", params.Pos->theta, params.Pos->Ntheta);
        for(int32_t i = 0; i < params.Pos->Ntheta; ++i) {
            params.Pos->theta[i] = RL(72.0) * (real)i;
        }
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        if constexpr(O3D) {
            ReadVector2(params, params.Pos->Ntheta, params.Pos->theta, ENVFile);
            CheckFix360Sweep(Pos->theta, Pos->Ntheta);
        } else {
            Default(params);
        }
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        ValidateVector2(
            params, params.Pos->Ntheta, params.Pos->theta, "Receiver bearings, theta");
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        Preprocess(params);
        EchoVector2(
            params, params.Pos->Ntheta, params.Pos->theta, RadDeg,
            "Receiver bearings, theta", "degrees");
    }
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const
    {
        // calculate angular spacing
        Pos->Delta_theta = FL(0.0);
        if(Pos->Ntheta != 1)
            Pos->Delta_theta = Pos->theta[Pos->Ntheta - 1] - Pos->theta[Pos->Ntheta - 2];

        trackallocate(params, "receiver bearing sin/cos table", Pos->t_rcvr, Pos->Ntheta);
        for(int32_t i = 0; i < params.Pos->Ntheta; ++i) {
            real theta            = DegRad * params.Pos->theta[i];
            params.Pos->t_rcvr[i] = vec2(STD::cos(theta), STD::sin(theta));
        }
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        trackdeallocate(params, params.Pos->theta);
        trackdeallocate(params, params.Pos->t_rcvr);
    }
};

} // namespace bhc
