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
 * Limits for tracing beams
 */
template<bool O3D> class BeamInfo : public ParamsModule<O3D> {
public:
    BeamInfo() {}
    virtual ~BeamInfo() {}

    virtual void SetupPre(bhcParams<O3D> &params) const override
    {
        BeamStructure<O3D> *Beam = params.Beam;
        memcpy(Beam->Type, "G S ", 4);
        MoveBeamType(params);

        Beam->rangeInKm  = true;
        Beam->autoDeltas = false;

        Beam->deltas = RL(0.0);
        Beam->Box.x = Beam->Box.y = RL(-1.0);
        if constexpr(O3D) Beam->Box.z = RL(-1.0);

        Beam->epsMultiplier = FL(1.0);
        Beam->rLoop         = RL(1.0);
        Beam->Nimage        = 1;
        Beam->iBeamWindow   = 4;
        Beam->Component     = 'P';
    }
    virtual void Default(bhcParams<O3D> &params) const override
    {
        BeamStructure<O3D> *Beam = params.Beam;
        if constexpr(O3D) {
            Beam->Box.x = Beam->Box.y = RL(101.0);
            Beam->Box.z               = RL(6000.0);
        } else {
            Beam->Box.y = RL(6000.0);
            Beam->Box.x = RL(101.0);
        }
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        BeamStructure<O3D> *Beam = params.Beam;

        LIST(ENVFile);
        ENVFile.Read(Beam->deltas);
        if constexpr(O3D) {
            ENVFile.Read(Beam->Box.x);
            ENVFile.Read(Beam->Box.y);
            ENVFile.Read(Beam->Box.z);
        } else {
            ENVFile.Read(Beam->Box.y); // LP: z then r
            ENVFile.Read(Beam->Box.x);
        }

        // *** Beam characteristics ***

        // no worry about the beam type if this is a ray trace run
        if(IsRayRun(Beam)) return;

        // Curvature change can cause overflow in grazing case
        // Suppress by setting BeamType[2] = 'Z'

        if(IsCervenyInfl(Beam)) {
            LIST(ENVFile);
            ENVFile.Read(&Beam->Type[1], 2);
            ENVFile.Read(Beam->epsMultiplier);
            ENVFile.Read(Beam->rLoop);

            // Images, windows
            LIST(ENVFile);
            ENVFile.Read(Beam->Nimage);
            ENVFile.Read(Beam->iBeamWindow);
            ENVFile.Read(Beam->Component);
        }
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        BeamStructure<O3D> *Beam = params.Beam;

        ENVFile << Beam->deltas;
        if constexpr(O3D) {
            ENVFile << (XYCOMP(Beam->Box) * RL(0.001)) << Beam->Box.z;
            ENVFile.write("! deltas (m), box X (km), box Y (km), box depth (m)\n");
        } else {
            ENVFile << Beam->Box.y << (Beam->Box.x * RL(0.001)); // LP: z then r
            ENVFile.write("! deltas (m), box depth (m), box range (km)\n");
        }

        if(IsRayRun(Beam)) return;

        if(IsCervenyInfl(Beam)) {
            ENVFile << std::string(&Beam->Type[1], 2);
            ENVFile << Beam->epsMultiplier << Beam->rLoop << (int32_t)0;
            ENVFile.write("! 'Width Curvature' epsMultiplier rLoop ISINGL (ignored)\n");

            ENVFile << Beam->Nimage << Beam->iBeamWindow
                    << std::string(&Beam->Component, 1);
            ENVFile.write("! Nimage iBeamWindow Component\n");
        }
    }
    virtual void Validate(bhcParams<O3D> &params) const override
    {
        BeamStructure<O3D> *Beam = params.Beam;
        Preprocess(params);

        bool boxerr = Beam->Box.x <= RL(0.0) || Beam->Box.y <= RL(0.0);
        if constexpr(O3D) boxerr = boxerr || Beam->Box.z <= RL(0.0);
        if(boxerr) { EXTERR("ReadEnvironment: Beam box not set up correctly"); }

        if(IsGeometricInfl(Beam) || IsSGBInfl(Beam)) {
            NULLSTATEMENT;
        } else if(IsCervenyInfl(Beam)) {
            switch(Beam->Type[2]) {
            case 'D': break;
            case 'Z': break;
            case 'S': break;
            default:
                EXTERR(
                    "ReadEnvironment: Unknown curvature condition '%c'", Beam->Type[2]);
            }
        } else {
            EXTERR(
                "ReadEnvironment: Unknown beam type (second letter of run type == '%c')",
                Beam->RunType[1]);
        }
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile    = GetInternal(params)->PRTFile;
        BeamStructure<O3D> *Beam = params.Beam;

        Preprocess(params);

        PRTFile
            << "\n_______________________________________________________________________"
               "___\n\n";
        PRTFile << std::setprecision(4);
        PRTFile << "\n Step length,       deltas = " << std::setw(11) << Beam->deltas
                << " m\n\n";
        if constexpr(O3D) {
            PRTFile << "Maximum ray x-range, Box.x  = " << std::setw(11) << Beam->Box.x
                    << " m\n";
            PRTFile << "Maximum ray y-range, Box.y  = " << std::setw(11) << Beam->Box.y
                    << " m\n";
            PRTFile << "Maximum ray z-range, Box.z  = " << std::setw(11) << Beam->Box.z
                    << " m\n";
        } else {
            PRTFile << "Maximum ray depth, Box.y  = " << std::setw(11) << Beam->Box.y
                    << " m\n";
            PRTFile << "Maximum ray range, Box.x  = " << std::setw(11)
                    << (Beam->Box.x * RL(0.001)) << "km\n";
        }

        // write info to prt file
        const char *tag = IsCervenyInfl(Beam) ? GetBeamWidthTag(Beam)
                                              : GetBeamTypeTag<O3D>(Beam);
        PRTFile << "\n" << tag << "\n";
        // LP: The first two aren't available here (only in run), and the third
        // is echoed below anyway.
        // PRTFile << "halfwidth  = " << halfwidth << "\n";
        // PRTFile << "epsilonOpt = " << epsilonOpt << "\n";
        // PRTFile << "EpsMult    = " << Beam->epsMultiplier << "\n\n";

        if(Beam->Type[3] == 'S') {
            PRTFile << "Beam shift in effect\n";
        } else {
            PRTFile << "No beam shift in effect\n";
        }

        if(!IsRayRun(Beam)) {
            if(IsCervenyInfl(Beam)) {
                PRTFile << "\n\nType of beam = " << Beam->Type[0] << "\n";
                switch(Beam->Type[2]) {
                case 'D': PRTFile << "Curvature doubling invoked\n"; break;
                case 'Z': PRTFile << "Curvature zeroing invoked\n"; break;
                case 'S': PRTFile << "Standard curvature condition\n"; break;
                }

                PRTFile << "Epsilon multiplier " << Beam->epsMultiplier << "\n";
                PRTFile << "Range for choosing beam width " << Beam->rLoop << "\n";

                // Images, windows
                PRTFile << "\nNumber of images, Nimage  = " << Beam->Nimage << "\n";
                PRTFile << "Beam windowing parameter  = " << Beam->iBeamWindow << "\n";
                PRTFile << "Component                 = " << Beam->Component << "\n";
            }

            PRTFile << "\n";
        }

        if constexpr(!O3D) {
            if(Beam->autoDeltas) {
                PRTFile << "\n Step length,       deltas = " << params.Beam->deltas
                        << " m (automatically selected)\n";
            }
        }
    }
    virtual void Preprocess(bhcParams<O3D> &params) const override
    {
        BeamStructure<O3D> *Beam = params.Beam;
        MoveBeamType(params);

        if(Beam->rangeInKm) {
            Beam->rangeInKm = false;
            Beam->Box.x *= FL(1000.0);                   // convert km to m
            if constexpr(O3D) Beam->Box.y *= FL(1000.0); // convert km to m
        }

        // Automatic step size selection
        if(Beam->deltas == FL(0.0)) {
            Beam->deltas = (params.Bdry->Bot.hs.Depth - params.Bdry->Top.hs.Depth)
                / FL(10.0);
            Beam->autoDeltas = true;
        }
    }

private:
    void MoveBeamType(bhcParams<O3D> &params) const
    {
        BeamStructure<O3D> *Beam = params.Beam;
        Beam->Type[3]            = Beam->RunType[6]; // selects beam shift option
        Beam->Type[0]            = Beam->RunType[1];
    }
};

}} // namespace bhc::module
