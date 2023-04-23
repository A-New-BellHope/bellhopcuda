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

template<bool O3D> class Title : public ParamsModule<O3D> {
public:
    Title() {}
    virtual ~Title() {}

    virtual void Default(bhcParams<O3D> &params) const override
    {
        SetTitle(params, "no env file");
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        std::string TempTitle;
        LIST(ENVFile);
        ENVFile.Read(TempTitle);
        SetTitle(params, TempTitle);
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        ENVFile << std::string(params.Title);
        ENVFile.write("! TITLE\n");
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        PRTFile << BHC_PROGRAMNAME "- " << params.Title << "\n";
    }

    inline void SetTitle(bhcParams<O3D> &params, const std::string &TempTitle) const
    {
        size_t l = bhc::min(sizeof(params.Title) - 1, TempTitle.size());
        memcpy(params.Title, TempTitle.c_str(), l);
        params.Title[l] = 0;
    }
};

}} // namespace bhc::module
