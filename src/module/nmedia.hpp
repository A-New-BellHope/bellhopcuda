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

template<bool O3D, bool R3D> class NMedia : public ParamsModule<O3D, R3D> {
public:
    NMedia() {}
    virtual ~NMedia() {}

    virtual void Default(bhcParams<O3D, R3D> &params) const override
    {
        params.atten->NMedia = 1;
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        LIST(ENVFile);
        ENVFile.Read(params.atten->NMedia);
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const override
    {
        if(params.atten->NMedia != 1) {
            EXTWARN("ReadEnvironment: Only one medium or layer is allowed in BELLHOP; "
                    "sediment layers must be handled using a reflection coefficient");
        }
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const override
    {
        PRTFile << "Dummy parameter NMedia = " << params.atten->NMedia << "\n";
    }
};

}} // namespace bhc::module