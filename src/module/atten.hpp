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

template<bool O3D, bool R3D> cpx crci(
    const bhcParams<O3D, R3D> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
extern template cpx crci<false, false>(
    const bhcParams<false, false> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
extern template cpx crci<true, false>(
    const bhcParams<true, false> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
extern template cpx crci<true, true>(
    const bhcParams<true, true> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);

template<bool O3D, bool R3D> class Atten : public ParamsModule<O3D, R3D> {
public:
    Atten() {}
    virtual ~Atten() {}

    virtual void SetupPre(bhcParams<O3D, R3D> &params) const
    {
        params.atten->t        = FL(20.0);
        params.atten->Salinity = FL(35.0);
        params.atten->pH       = FL(8.0);
        params.atten->z_bar    = FL(0.0);
    }
};

}} // namespace bhc::module
