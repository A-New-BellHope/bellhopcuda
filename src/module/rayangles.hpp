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

template<bool O3D, bool BEARING> class RayAngles : public ParamsModule<O3D> {
public:
    RayAngles() {}
    virtual ~RayAngles() {}

    virtual void Init(bhcParams<O3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        a.angles     = nullptr;
    }

    virtual void SetupPre(bhcParams<O3D> &params) const override
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
        a.inDegrees = true;
    }

    virtual void Default(bhcParams<O3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        if constexpr(BEARING) {
            if constexpr(O3D) { a.n = 5; }
            trackallocate(params, FuncName, a.angles, a.n);
            for(int32_t i = 0; i < a.n; ++i) {
                a.angles[i] = (float)(i * 360) / (float)(a.n);
            }
        } else {
            EstimateNumAngles(params);
            trackallocate(params, FuncName, a.angles, a.n);
            if(a.n < 3) EXTERR("Internal error in default RayAnglesElevation setup");
            a.angles[0] = RL(-20.0);
            a.angles[1] = RL(20.0);
            a.angles[2] = FL(-999.9);
            SubTab(a.angles, a.n);
        }
    }

    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        if constexpr(BEARING && !O3D) {
            Default(params);
        } else {
            AngleInfo &a = GetAngle(params);
            LIST(ENVFile);
            ENVFile.Read(a.n);
            if(params.Bdry->Top.hs.Opt[5] == 'I') {
                // option to trace a single beam
                ENVFile.Read(a.iSingle);
            }
            if(a.n == 0) EstimateNumAngles(params);
            trackallocate(params, FuncName, a.angles, bhc::max(3, a.n));

            if(a.n > 2) a.angles[2] = FL(-999.9);
            LIST(ENVFile);
            ENVFile.Read(a.angles, a.n);
            SubTab(a.angles, a.n);
            Sort(a.angles, a.n);
            [[maybe_unused]] bool dummy;
            CheckFix360Sweep(a.angles, a.n, dummy);
        }
    }

    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        if constexpr(!BEARING || O3D) {
            AngleInfo &a = GetAngle(params);
            const char *nLabel, *xLabel;
            if constexpr(BEARING) {
                nLabel = "Nbeta";
                xLabel = "beta (degrees) bearing angle fan";
            } else {
                nLabel = "NBEAMS";
                xLabel = "ALPHA (degrees)";
            }
            UnSubTab(ENVFile, a.angles, a.n, nLabel, xLabel, RadDeg);
        }
    }

    void ExtSetup(bhcParams<O3D> &params, int32_t n) const
    {
        AngleInfo &a = GetAngle(params);
        a.n          = n;
        trackallocate(params, FuncName, a.angles, a.n);
    }

    virtual void Validate(bhcParams<O3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        if constexpr(BEARING && !O3D) {
            if(a.n != 1 || a.angles[0] != RL(0.0)) {
                EXTERR("Invalid beam bearing angles setup for 2D");
            }
        } else {
            ValidateVector(params, a.angles, a.n, FuncName);

            if(a.n > 1 && a.angles[a.n - 1] == a.angles[0]) {
                EXTERR("%s: First and last beam take-off angle are identical", FuncName);
            }
            if(params.Bdry->Top.hs.Opt[5] == 'I' && (a.iSingle < 1 || a.iSingle > a.n)) {
                EXTERR("%s: Selected beam, iSingle not in [1, a.n]", FuncName);
            }
        }
    }

    virtual void Echo(bhcParams<O3D> &params) const override
    {
        if constexpr(!BEARING || O3D) {
            PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
            AngleInfo &a          = GetAngle(params);
            Preprocess(params);
            if constexpr(!BEARING) {
                PRTFile << "_____________________________________________________________"
                           "________"
                           "_____\n";
            }
            PRTFile << "\n   Number of beams in " << (BEARING ? "bearing  " : "elevation")
                    << "   = " << a.n << "\n";
            if(a.iSingle > 0) PRTFile << "Trace only beam number " << a.iSingle << "\n";
            PRTFile << "   Beam take-off angles (degrees)\n";
            EchoVector(a.angles, a.n, PRTFile, 10, "", RadDeg);
        }
    }

    virtual void Preprocess(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        AngleInfo &a          = GetAngle(params);
        if(a.n <= 0) { EXTERR("%s number of beams = %d is invalid", FuncName, a.n); }

        if constexpr(BEARING) {
            if(GetInternal(params)->dim == 4 && !IsRayRun(params.Beam)) {
                // Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
                PRTFile
                    << "\nReplacing beam take-off angles, beta, with receiver bearing "
                       "lines, theta\n";
                trackdeallocate(params, a.angles);

                a.n         = params.Pos->Ntheta;
                a.inDegrees = true;
                trackallocate(
                    params, "beam angles replaced with receiver bearing angles", a.angles,
                    a.n);
                for(int32_t i = 0; i < a.n; ++i)
                    a.angles[i] = params.Pos->theta[i]; // a.n should = params.Pos->Ntheta
            }
        }

        if(a.inDegrees) {
            // convert to radians
            for(int32_t i = 0; i < a.n; ++i) a.angles[i] *= DegRad;
            a.inDegrees = false;
        }

        // angular spacing between beams
        a.d = FL(0.0);
        if(a.n != 1) a.d = (a.angles[a.n - 1] - a.angles[0]) / (a.n - 1);
    }

    virtual void Finalize(bhcParams<O3D> &params) const override
    {
        AngleInfo &a = GetAngle(params);
        trackdeallocate(params, a.angles);
    }

private:
    constexpr static const char *FuncName = BEARING ? "RayAnglesBearing"
                                                    : "RayAnglesElevation";

    inline AngleInfo &GetAngle(bhcParams<O3D> &params) const
    {
        if constexpr(BEARING)
            return params.Angles->beta;
        else
            return params.Angles->alpha;
    }
    /**
     * automatically estimate n to use
     */
    inline void EstimateNumAngles(bhcParams<O3D> &params) const
    {
        AngleInfo &a = GetAngle(params);
        if(IsRayRun(params.Beam)) {
            a.n = 50; // For a ray trace plot, we don't want too many rays ...
        } else {
            // you're letting ME choose? OK: ideas based on an isospeed ocean
            // limit based on phase of adjacent beams at maximum range
            constexpr real c0 = FL(1500.0);
            a.n = bhc::max(
                (int)((BEARING ? FL(0.1) : FL(0.3)) * params.Pos->Rr[params.Pos->NRr - 1] 
                * params.freqinfo->freq0 / c0),
                300);

            if constexpr(!BEARING) {
                // limit based on having beams that are thin with respect to the water
                // depth assumes also a full 360 degree angular spread of rays
                // Should check which Depth is used here, in case where there is a
                // variable bathymetry
                real Depth = params.Bdry->Bot.hs.Depth - params.Bdry->Top.hs.Depth;
                real d_theta_recommended = STD::atan(
                    Depth / (FL(10.0) * params.Pos->Rr[params.Pos->NRr - 1]));
                a.n = bhc::max((int)(REAL_PI / d_theta_recommended), a.n);
            }
        }
    }
};

template<bool O3D> using RayAnglesElevation = RayAngles<O3D, false>;
template<bool O3D> using RayAnglesBearing   = RayAngles<O3D, true>;

}} // namespace bhc::module
