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

namespace bhc { namespace mode {

/**
 * Like ParamsModule, but for outputs, and fewer steps.
 */
template<bool O3D, bool R3D> class ModeModule {
public:
    ModeModule() {}
    virtual ~ModeModule() {}

    /// Initialization and defaults.
    virtual void Init(bhcOutputs<O3D, R3D> &) const {}
    /// Preprocessing as part of run.
    virtual void Preprocess(bhcParams<O3D> &, bhcOutputs<O3D, R3D> &) const {}
    /// Run the simulation.
    virtual void Run(bhcParams<O3D> &, bhcOutputs<O3D, R3D> &) const = 0;
    /// Postprocess after run is complete.
    virtual void Postprocess(bhcParams<O3D> &, bhcOutputs<O3D, R3D> &) const {}
    /// Write results to disk.
    virtual void Writeout(const bhcParams<O3D> &, const bhcOutputs<O3D, R3D> &) const {}
    /// Read results from disk.
    virtual void Readout(bhcParams<O3D> &, bhcOutputs<O3D, R3D> &, const char *) const {}
    /// Deallocate memory.
    virtual void Finalize(bhcParams<O3D> &, bhcOutputs<O3D, R3D> &) const {}
};

}} // namespace bhc::mode
