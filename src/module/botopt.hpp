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
 *
 */
template<bool O3D, bool R3D> class BotOpt {
public:
    BotOpt() {}
    virtual ~BotOpt() {}

    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        memcpy(params.Bdry->Bot.hs.Opt, "R-    ", 6);
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.Bdry->Bot.hs.Opt, 6); // LP: LDIFile fills rest with ' '
        ENVFile.Read(params.Bdry->Bot.hsx.Sigma);
    }
    virtual void SetupPost(bhcParams<O3D, R3D> &params) const override
    {
        params.Bdry->Bot.hs.bc = params.Bdry->Bot.hs.Opt[0];
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

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
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
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
