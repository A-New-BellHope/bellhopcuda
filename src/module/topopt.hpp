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
#include "boundarycond.hpp"

namespace bhc { namespace module {

/**
 * LP: Read top halfspace options; 4 out of the 6 entries are general program
 * options.
 */
template<bool O3D> class TopOpt : public ParamsModule<O3D> {
public:
    TopOpt() {}
    virtual ~TopOpt() {}

    virtual void Default(bhcParams<O3D> &params) const override
    {
        // SSP (clinear), top bc (vacuum), atten units (dB/wavelength),
        // add vol atten (none), altimetry (none), dev mode (off)
        memcpy(params.Bdry->Top.hs.Opt, "CVW - ", 6);
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.Bdry->Top.hs.Opt, 6); // LP: LDIFile fills rest with ' '

        SetupPost(params);

        // optional addition of volume attenuation using standard formulas

        AttenInfo *atten = params.atten;
        switch(params.ssp->AttenUnit[1]) {
        case 'F':
            LIST(ENVFile);
            ENVFile.Read(atten->t);
            ENVFile.Read(atten->Salinity);
            ENVFile.Read(atten->pH);
            ENVFile.Read(atten->z_bar);
            break;
        case 'B':
            LIST(ENVFile);
            ENVFile.Read(atten->NBioLayers);
            for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio) {
                LIST(ENVFile);
                ENVFile.Read(atten->bio[iBio].z1);
                ENVFile.Read(atten->bio[iBio].z2);
                ENVFile.Read(atten->bio[iBio].f0);
                ENVFile.Read(atten->bio[iBio].q);
                ENVFile.Read(atten->bio[iBio].a0);
            }
        }
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        ENVFile << std::string(params.Bdry->Top.hs.Opt, 6);
        ENVFile.write("! SSP (");
        switch(params.Bdry->Top.hs.Opt[0]) {
        case 'N': ENVFile.write("n2linear"); break;
        case 'C': ENVFile.write("clinear"); break;
        case 'P': ENVFile.write("pchip"); break;
        case 'S': ENVFile.write("spline"); break;
        case 'Q': ENVFile.write("quad"); break;
        case 'H': ENVFile.write("hexahedral"); break;
        case 'A': ENVFile.write("analytic"); break;
        default: ENVFile.write("error!");
        }
        ENVFile.write("), top bc (");
        BoundaryCond<O3D, true>::WriteBCTag(params.Bdry->Top.hs.Opt[1], ENVFile);
        ENVFile.write("), atten units, add vol atten, altimetry, dev mode / broadband\n");
    }
    virtual void SetupPost(bhcParams<O3D> &params) const override
    {
        params.ssp->Type         = params.Bdry->Top.hs.Opt[0];
        params.Bdry->Top.hs.bc   = params.Bdry->Top.hs.Opt[1];
        params.ssp->AttenUnit[0] = params.Bdry->Top.hs.Opt[2];
        params.ssp->AttenUnit[1] = params.Bdry->Top.hs.Opt[3];
    }
    virtual void Validate(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

        // Copy values back to TopOpt in case the user modified the other ones
        params.Bdry->Top.hs.Opt[0] = params.ssp->Type;
        params.Bdry->Top.hs.Opt[1] = params.Bdry->Top.hs.bc;
        params.Bdry->Top.hs.Opt[2] = params.ssp->AttenUnit[0];
        params.Bdry->Top.hs.Opt[3] = params.ssp->AttenUnit[1];

        // SSP approximation options

        switch(params.ssp->Type) {
        case 'N': break;
        case 'C': break;
        case 'P': break;
        case 'S': break;
        case 'Q': {
            // LP: This just checks for existence, moved actual open for reading
            // to InitQuad.
            std::ifstream SSPFile;
            SSPFile.open(GetInternal(params)->FileRoot + ".ssp");
            if(!SSPFile.good()) {
                PRTFile << "SSPFile = " << GetInternal(params)->FileRoot << ".ssp\n";
                EXTERR(BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file");
            }
        } break;
        case 'H': {
            // LP: This just checks for existence, moved actual open for reading
            // to InitHexahedral.
            std::ifstream SSPFile;
            SSPFile.open(GetInternal(params)->FileRoot + ".ssp");
            if(!SSPFile.good()) {
                PRTFile << "SSPFile = " << GetInternal(params)->FileRoot << ".ssp\n";
                EXTERR(BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file");
            }
        } break;
        case 'A': break;
        default: EXTERR("ReadEnvironment: Unknown option for SSP approximation");
        }

        // Attenuation options

        switch(params.ssp->AttenUnit[0]) {
        case 'N': break;
        case 'F': break;
        case 'M': break;
        case 'W': break;
        case 'Q': break;
        case 'L': break;
        default: EXTERR("ReadEnvironment: Unknown attenuation units");
        }

        // optional addition of volume attenuation using standard formulas

        switch(params.ssp->AttenUnit[1]) {
        case 'T': break;
        case 'F': break;
        case 'B': break;
        case ' ': break;
        default: EXTERR("ReadEnvironment: Unknown top option letter in fourth position");
        }

        switch(params.Bdry->Top.hs.Opt[4]) {
        case '~':
        case '*': break;
        case '-':
        case '_':
        case ' ': break;
        default: EXTERR("ReadEnvironment: Unknown top option letter in fifth position");
        }

        switch(params.Bdry->Top.hs.Opt[5]) {
        case 'I': break;
        case ' ': break;
        case 'B':
            EXTERR("ReadEnvironment: BELLHOP/BELLHOP3D does not properly support "
                   "wideband runs / freqvec / TopOpt[5] == 'B'");
            [[fallthrough]]; // LP: break here would be unreachable
        default: EXTERR("ReadEnvironment: Unknown top option letter in sixth position\n");
        }
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        PRTFile << "\n";

        // SSP approximation options

        switch(params.ssp->Type) {
        case 'N': PRTFile << "    N2-linear approximation to SSP\n"; break;
        case 'C': PRTFile << "    C-linear approximation to SSP\n"; break;
        case 'P': PRTFile << "    PCHIP approximation to SSP\n"; break;
        case 'S': PRTFile << "    Spline approximation to SSP\n"; break;
        case 'Q': PRTFile << "    Quad approximation to SSP\n"; break;
        case 'H': PRTFile << "    Hexahedral approximation to SSP\n"; break;
        case 'A': PRTFile << "    Analytic SSP option\n"; break;
        }

        // Attenuation options

        switch(params.ssp->AttenUnit[0]) {
        case 'N': PRTFile << "    Attenuation units: nepers/m\n"; break;
        case 'F': PRTFile << "    Attenuation units: dB/mkHz\n"; break;
        case 'M': PRTFile << "    Attenuation units: dB/m\n"; break;
        case 'W': PRTFile << "    Attenuation units: dB/wavelength\n"; break;
        case 'Q': PRTFile << "    Attenuation units: Q\n"; break;
        case 'L': PRTFile << "    Attenuation units: Loss parameter\n"; break;
        }

        // optional addition of volume attenuation using standard formulas

        AttenInfo *atten = params.atten;
        switch(params.ssp->AttenUnit[1]) {
        case 'T': PRTFile << "    THORP volume attenuation added\n"; break;
        case 'F':
            PRTFile << "    Francois-Garrison volume attenuation added\n";
            PRTFile << std::setprecision(4);
            PRTFile << " T = " << std::setw(11) << atten->t
                    << "degrees   S = " << std::setw(11) << atten->Salinity
                    << " psu   pH = " << std::setw(11) << atten->pH
                    << " z_bar = " << std::setw(11) << " m\n";
            break;
        case 'B':
            PRTFile << "    Biological attenuation\n";
            PRTFile << "      Number of Bio Layers = " << atten->NBioLayers << "\n";
            for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio) {
                PRTFile << "      Top    of layer = " << atten->bio[iBio].z1 << " m\n";
                PRTFile << "      Bottom of layer = " << atten->bio[iBio].z2 << " m\n";
                PRTFile << "      Resonance frequency = " << atten->bio[iBio].f0
                        << " Hz\n";
                PRTFile << "      Q = " << atten->bio[iBio].q << "\n";
                PRTFile << "      a0 = " << atten->bio[iBio].a0 << "\n";
            }
            break;
        }

        switch(params.Bdry->Top.hs.Opt[4]) {
        case '~':
        case '*': PRTFile << "    Altimetry file selected\n"; break;
        }

        if(params.Bdry->Top.hs.Opt[5] == 'I') {
            PRTFile << "    Development options enabled\n";
        }
    }
};

}} // namespace bhc::module
