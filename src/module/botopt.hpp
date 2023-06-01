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

template<bool O3D> class BotOpt : public ParamsModule<O3D> {
public:
    BotOpt() {}
    virtual ~BotOpt() {}

    virtual void Default(bhcParams<O3D> &params) const override
    {
        memcpy(params.Bdry->Bot.hs.Opt, "R-    ", 6);
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.Bdry->Bot.hs.Opt, 6); // LP: LDIFile fills rest with ' '
        ENVFile.Read(params.Bdry->Bot.hsx.Sigma);
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        ENVFile << std::string(params.Bdry->Bot.hs.Opt, 6) << params.Bdry->Bot.hsx.Sigma;
        ENVFile.write("! bot bc (");
        BoundaryCond<O3D, false>::WriteBCTag(params.Bdry->Bot.hs.Opt[0], ENVFile);
        ENVFile.write("), bathymetry, 4 spaces; Sigma (printed but ignored)\n");
    }
    virtual void SetupPost(bhcParams<O3D> &params) const override
    {
        params.Bdry->Bot.hs.bc = params.Bdry->Bot.hs.Opt[0];
    }
    virtual void Validate(bhcParams<O3D> &params) const override
    {
        switch(params.Bdry->Bot.hs.Opt[1]) {
        case '~':
        case '*': break;
        case '-':
        case '_':
        case ' ': break;
        default:
            EXTERR(
                "Unknown bottom option letter in second position: Bdry->Bot.hs.Opt[1] == "
                "'%c'",
                params.Bdry->Bot.hs.Opt[1]);
        }
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

        PRTFile << "\n RMS roughness = " << std::setw(10) << std::setprecision(3)
                << params.Bdry->Bot.hsx.Sigma << "\n";

        switch(params.Bdry->Bot.hs.Opt[1]) {
        case '~':
        case '*': PRTFile << "    Bathymetry file selected\n";
        }
    }
};

}} // namespace bhc::module
