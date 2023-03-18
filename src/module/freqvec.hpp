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
#include "paramsmodule.hpp"

namespace bhc {

/**
 * Optionally reads a vector of source frequencies for a broadband run
 * If the broadband option is not selected, then the input freq (a scalar) is stored in
 * the frequency vector
 */
template<bool O3D, bool R3D> class FreqVec : public ParamsModule {
public:
    FreqVec() {}
    virtual ~FreqVec() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        params.freqinfo->freqVec = nullptr;
    }
    virtual void SetupPre(bhcParams<O3D, R3D> &params) const
    {
        params.freqinfo->Nfreq = 1;
    }
    virtual void Default(bhcParams<O3D, R3D> &params) const
    {
        trackallocate(
            params, "default source frequencies", freqinfo->freqVec,
            params.freqinfo->Nfreq);
        freqinfo->freqVec[0] = freqinfo->freq0;
    }
    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        if(params.Bdry->Top.hs.Opt[5] == 'B') {
            ReadVector2(
                params, params.freqinfo->Nfreq, params.freqinfo->freqVec, ENVFile);
        } else {
            Default(params);
        }
    }
    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        ValidateVector2(
            params, params.freqinfo->Nfreq, params.freqinfo->freqVec, "frequencies");
    }
    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        Preprocess(params);
        EchoVector2(
            params, params.freqinfo->Nfreq, params.freqinfo->freqVec, RL(1.0),
            "Frequencies", "Hz");
    }
    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        trackdeallocate(params, params.freqinfo->freqVec);
    }
};

} // namespace bhc
