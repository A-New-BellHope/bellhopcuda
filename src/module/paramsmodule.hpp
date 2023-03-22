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

namespace bhc { namespace module {

/**
 * Child classes are responsible for the initialization, defaults, reading,
 * etc. of some set of parameters within params. In addition to the methods
 * here, modules may provide a ExtSetup method to, for example, allocate an
 * array of a given size from the external API. Modules should not contain any
 * member variables--this is effectively a "static virtual" class, but this
 * concept is not supported in C++.
 */
template<bool O3D, bool R3D> class ParamsModule {
public:
    virtual ~ParamsModule() {}

    /// True one-time initialization, e.g. set pointers to arrays to nullptr.
    virtual void Init(bhcParams<O3D, R3D> &params) const {}
    /// Called before Default or Read for common setup. Can set some defaults
    /// here which are needed in case the env file does not write certain
    /// variables.
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const {}
    /// Set the parameters to some reasonable default values in place of Read.
    virtual void Default(bhcParams<O3D, R3D> &params) const {}
    /// Read the parameters from the environment file or other input files.
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {}
    /// Called after Default or Read for common setup.
    virtual void SetupPost(bhcParams<O3D, R3D> &params) const {}
    /// Check if the parameters are valid values. Throws errors if not.
    virtual void Validate(const bhcParams<O3D, R3D> &params) const {}
    /// Writes info about the parameters to the print file emulator.
    virtual void Echo(const bhcParams<O3D, R3D> &params) const {}
    /// Modifies the parameters before processing, e.g. km to m. Module must add
    /// flags to params to track whether this has been done or not.
    virtual void Preprocess(bhcParams<O3D, R3D> &params) const {}
    /// Deallocate memory.
    virtual void Finalize(bhcParams<O3D, R3D> &params) const {}

private:
    ParamsModule() {}
};

}} // namespace bhc::module
