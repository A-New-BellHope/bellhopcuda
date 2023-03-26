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

template<bool O3D, bool R3D> class Title : public ParamsModule<O3D, R3D> {
public:
    Title() {}
    virtual ~Title() {}

    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        SetTitle(params, "no env file");
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        std::string TempTitle;
        LIST(ENVFile);
        ENVFile.Read(TempTitle);
        SetTitle(params, TempTitle);
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        PRTFile << params.Title << "\n";
    }

private:
    inline void SetTitle(bhcParams<O3D, R3D> &params, const std::string &TempTitle) const
    {
        // Prepend model name to title
        TempTitle = BHC_PROGRAMNAME "- " + TempTitle;
        size_t l  = bhc::min(sizeof(params.Title) - 1, TempTitle.size());
        memcpy(params.Title, TempTitle.c_str(), l);
        params.Title[l] = 0;
    }
};

}} // namespace bhc::module
