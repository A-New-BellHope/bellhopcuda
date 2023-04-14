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

template<bool O3D> cpx crci(
    const bhcParams<O3D> &params, real z, real c, real alpha, const char (&AttenUnit)[2]);
extern template cpx crci<false>(
    const bhcParams<false> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
extern template cpx crci<true>(
    const bhcParams<true> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);

template<bool O3D> class Atten : public ParamsModule<O3D> {
public:
    Atten() {}
    virtual ~Atten() {}

    virtual void SetupPre(bhcParams<O3D> &params) const override
    {
        params.atten->t        = FL(20.0);
        params.atten->Salinity = FL(35.0);
        params.atten->pH       = FL(8.0);
        params.atten->z_bar    = FL(0.0);
    }
    virtual void Default(bhcParams<O3D> &) const override {}
};

}} // namespace bhc::module
