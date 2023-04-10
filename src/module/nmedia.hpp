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

template<bool O3D> class NMedia : public ParamsModule<O3D> {
public:
    NMedia() {}
    virtual ~NMedia() {}

    virtual void Default(bhcParams<O3D> &params) const override
    {
        params.atten->NMedia = 1;
    }
    virtual void Read(bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.atten->NMedia);
    }
    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        ENVFile << params.atten->NMedia;
        ENVFile.write("! NMEDIA\n");
    }
    virtual void Validate(bhcParams<O3D> &params) const override
    {
        if(params.atten->NMedia != 1) {
            EXTWARN("ReadEnvironment: Only one medium or layer is allowed in BELLHOP; "
                    "sediment layers must be handled using a reflection coefficient");
        }
    }
    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        PRTFile << "Dummy parameter NMedia = " << params.atten->NMedia << "\n";
    }
};

}} // namespace bhc::module
