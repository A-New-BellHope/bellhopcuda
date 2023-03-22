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
#include "misc.hpp"

namespace bhc { namespace mode {

template<bool O3D, bool R3D> class Eigen : public Field {
public:
    Eigen() {}
    virtual ~Eigen() {}

    virtual void Init(bhcOutputs<O3D, R3D> &outputs) const
    {
        outputs.eigen->hits    = nullptr;
        outputs.eigen->neigen  = 0;
        outputs.eigen->memsize = 0;
    }

    virtual void Preprocess(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {
        Field::Preprocess(params);
        TODO();
    }

    virtual void Postprocess(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {}

    virtual void Writeout(
        bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs) const
    {}

    virtual void Finalize(bhcOutputs<O3D, R3D> &outputs) const
    {
        trackdeallocate(params, outputs.eigen->hits);
    }
};

}} // namespace bhc::mode
