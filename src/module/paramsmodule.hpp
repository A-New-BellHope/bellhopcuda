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

namespace bhc { namespace module {

/**
 * Child classes are responsible for the initialization, defaults, reading,
 * etc. of some set of parameters within params. In addition to the methods
 * here, modules may provide a ExtSetup method to, for example, allocate an
 * array of a given size from the external API. Modules should not contain any
 * member variables--this is effectively a "static virtual" class, but this
 * concept is not supported in C++.
 */
template<bool O3D> class ParamsModule {
public:
    ParamsModule() {}
    virtual ~ParamsModule() {}

    /// True one-time initialization, e.g. set pointers to arrays to nullptr.
    virtual void Init(bhcParams<O3D> &) const {}
    /// Called before Default or Read for common setup. Can set some defaults
    /// here which are needed in case the env file does not write certain
    /// variables.
    virtual void SetupPre(bhcParams<O3D> &) const {}
    /// Set the parameters to some reasonable default values in place of Read.
    virtual void Default(bhcParams<O3D> &) const = 0;
    /// Read the parameters from the environment file or other input files.
    virtual void Read(bhcParams<O3D> &, LDIFile &, HSInfo &) const {}
    /// Write the parameters to an environment file and/or other output files.
    virtual void Write(bhcParams<O3D> &, LDOFile &) const {}
    /// Called after Default or Read for common setup.
    virtual void SetupPost(bhcParams<O3D> &) const {}
    /// Check if the parameters are valid values. Throws errors if not.
    virtual void Validate(bhcParams<O3D> &) const {}
    /// Writes info about the parameters to the print file emulator.
    virtual void Echo(bhcParams<O3D> &) const {}
    /// Modifies the parameters before processing, e.g. km to m. Module must add
    /// flags to params to track whether this has been done or not.
    virtual void Preprocess(bhcParams<O3D> &) const {}
    /// Deallocate memory.
    virtual void Finalize(bhcParams<O3D> &) const {}
};

}} // namespace bhc::module
