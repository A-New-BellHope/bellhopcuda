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

template<bool O3D, bool R3D, bool BEARING> class RayAngles : public ParamsModule {
public:
    RayAngles() {}
    virtual ~RayAngles() {}

    const char *GetFuncName() const
    {
        if constexpr(BEARING)
            return "RayAnglesBearing";
        else
            return "RayAnglesElevation";
    }
    AngleInfo &GetAngle(bhcParams<O3D, R3D> &params) const
    {
        if constexpr(BEARING)
            return params.Angles->beta;
        else
            return params.Angles->alpha;
    }
    /**
     * automatically estimate n to use
     */
    void EstimateNumAngles(bhcParams<O3D, R3D> &params, AngleInfo &a) const
    {
        if(IsRayRun(params.Beam)) {
            a.n = 50; // For a ray trace plot, we don't want too many rays ...
        } else {
            // you're letting ME choose? OK: ideas based on an isospeed ocean
            // limit based on phase of adjacent beams at maximum range
            a.n = bhc::max(
                (int)((BEARING ? FL(0.1) : FL(0.3)) * params.Pos->Rr[params.Pos->NRr - 1] 
                * params.freqinfo->freq0 / c0),
                300);

            if constexpr(!BEARING) {
                // limit based on having beams that are thin with respect to the water
                // depth assumes also a full 360 degree angular spread of rays Should
                // check which Depth is used here, in case where there is a variable
                // bathymetry
                real d_theta_recommended = STD::atan(
                    Depth / (FL(10.0) * params.Pos->Rr[params.Pos->NRr - 1]));
                a.n = bhc::max((int)(REAL_PI / d_theta_recommended), a.n);
            }
        }
    }

    virtual void Init(bhcParams<O3D, R3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        a.angles     = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        if constexpr(BEARING) {
            a.n = 1; // LP: Don't pick automatically, just 1 angle
        } else {
            a.n = 0; // LP: Pick automatically
        }
        // LP: not a typo; this is an index, one less than the start of the array,
        // which in Fortran (and in the env file!) is 0. This gets converted to 0-
        // indexed when it is used.
        a.iSingle   = 0;
        a.InDegrees = true;
    }
    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        if constexpr(BEARING) {
            if constexpr(O3D) { a.n = 5; }
            trackallocate(params, GetFuncName(), a.angles, a.n);
            for(int32_t i = 0; i < a.n; ++i) {
                a.angles[i] = (float)(i * 360) / (float)(a.n);
            }
        } else {
            EstimateNumAngles(params, a);
            trackallocate(params, GetFuncName(), a.angles, a.n);
            if(a.n < 3) EXTERR("Internal error in default RayAnglesElevation setup");
            a.angles[0] = RL(-20.0);
            a.angles[1] = RL(20.0);
            a.angles[2] = FL(-999.9);
            SubTab(a.angles, a.n);
        }
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        AngleInfo &a = GetAngle(params);

        LIST(ENVFile);
        ENVFile.Read(a.n);
        if(params.Bdry->Top.hs.Opt[5] == 'I') {
            // option to trace a single beam
            ENVFile.Read(a.iSingle);
        }
        if(a.n == 0) { EstimateNumAngles(params, a); }
        trackallocate(params, GetFuncName(), a.angles, bhc::max(3, a.n));

        if(a.n > 2) a.angles[2] = FL(-999.9);
        LIST(ENVFile);
        ENVFile.Read(a.angles, a.n);
        SubTab(a.angles, a.n);
        Sort(a.angles, a.n);
        CheckFix360Sweep(a.angles, a.n);
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        ValidateVector2(params, a.n, a.angles, GetFuncName());

        if(a.n > 1 && a.angles[a.n - 1] == a.angles[0]) {
            EXTERR("%s: First and last beam take-off angle are identical", GetFuncName());
        }
        if(params.Bdry->Top.hs.Opt[5] == 'I' && (a.iSingle < 1 || a.iSingle > a.n)) {
            EXTERR("%s: Selected beam, iSingle not in [1, a.n]", GetFuncName());
        }
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        Preprocess(params);
        if constexpr(!BEARING) {
            PRTFile
                << "_____________________________________________________________________"
                   "_____\n";
        }
        PRTFile << "\n   Number of beams in " << (BEARING ? "bearing  " : "elevation")
                << "   = " << a.n << "\n";
        if(a.iSingle > 0) PRTFile << "Trace only beam number " << a.iSingle << "\n";
        PRTFile << "   Beam take-off angles (degrees)\n";
        EchoVector(a.angles, a.n, RadDeg, PRTFile);
    }
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        if constexpr(BEARING && O3D && !R3D) {
            // Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
            if(!IsRayRun(params.Beam)) {
                PRTFile
                    << "\nReplacing beam take-off angles, beta, with receiver bearing "
                       "lines, theta\n";
                trackdeallocate(params, a.angles);

                a.n = params.Pos->Ntheta;
                trackallocate(
                    params, "beam angles replaced with receiver bearing angles", a.angles,
                    a.n);
                for(int32_t i = 0; i < a.n; ++i)
                    a.angles[i] = params.Pos->theta[i]; // a.n should = params.Pos->Ntheta
            }
        }

        if(a.InDegrees) {
            // convert to radians
            for(int32_t i = 0; i < a.n; ++i) a.angles[i] *= DegRad;
            a.InDegrees = false;
        }

        // angular spacing between beams
        a.d = FL(0.0);
        if(a.n != 1) a.d = (a.angles[a.n - 1] - a.angles[0]) / (a.n - 1);
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        trackdeallocate(params, a.angles);
    }
};

template<bool O3D, bool R3D> using RayAnglesElevation = RayAngles<O3D, R3D, false>;
template<bool O3D, bool R3D> using RayAnglesBearing   = RayAngles<O3D, R3D, true>;

}} // namespace bhc::module
