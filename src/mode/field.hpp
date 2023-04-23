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
#include "modemodule.hpp"
#include "fieldimpl.hpp"
#include "../influence.hpp"

namespace bhc { namespace mode {

/**
 * Parent class for field modes (TL, eigen, arr).
 */
template<bool O3D, bool R3D> class Field : public ModeModule<O3D, R3D> {
public:
    Field() {}
    virtual ~Field() {}

    virtual void Preprocess(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &) const override
    {
        PreRun_Influence<O3D, R3D>(params);
    }

    virtual void Run(bhcParams<O3D> &params, bhcOutputs<O3D, R3D> &outputs) const override
    {
        RunFieldModesSelInfl<O3D, R3D>(params, outputs);
    }
};

}} // namespace bhc::mode
