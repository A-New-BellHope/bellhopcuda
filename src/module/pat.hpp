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
 * Source Beam Pattern
 */
template<bool O3D, bool R3D> class Pat {
public:
    Pat() {}
    virtual ~Pat() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        params.beaminfo->SrcBmPat = nullptr;
    }

    virtual void SetupPre(bhcParams<O3D, R3D> &params) const
    {
        params.beaminfo->SBPFlag = params.Beam->RunType[2];
        params.beaminfo->SBPIndB = true;
    }

    virtual void Default(bhcParams<O3D, R3D> &params) const
    {
        BeamInfo *beaminfo = params.beaminfo;
        beaminfo->NSBPPts  = 2;
        trackallocate(
            params, "default/trivial source beam pattern", beaminfo->SrcBmPat, 2 * 2);
        beaminfo->SrcBmPat[0 * 2 + 0] = FL(-180.0);
        beaminfo->SrcBmPat[0 * 2 + 1] = FL(0.0);
        beaminfo->SrcBmPat[1 * 2 + 0] = FL(180.0);
        beaminfo->SrcBmPat[1 * 2 + 1] = FL(0.0);
    }

    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        BeamInfo *beaminfo = params.beaminfo;
        if(beaminfo->SBPFlag != '*') {
            Default(params);
            return;
        }

        LDIFile SBPFile(GetInternal(params), GetInternal(params)->FileRoot, ".sbp");
        if(!SBPFile.Good()) {
            PRTFile << "SBPFile = " << GetInternal(params)->FileRoot << ".sbp\n";
            EXTERR("BELLHOP-ReadPat: Unable to open source beampattern file");
        }

        LIST(SBPFile);
        SBPFile.Read(beaminfo->NSBPPts);
        trackallocate(
            params, "source beam pattern", beaminfo->SrcBmPat, beaminfo->NSBPPts * 2);

        for(int32_t i = 0; i < beaminfo->NSBPPts; ++i) {
            LIST(SBPFile);
            SBPFile.Read(&beaminfo->SrcBmPat[2 * i], 2);
        }
    }

    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        BeamInfo *beaminfo = params.beaminfo;
        if(!monotonic(beaminfo->SrcBmPat, beaminfo->NSBPPts, 2, 0)) {
            EXTERR("BELLHOP-ReadPat: Source beam pattern angles are not monotonic");
        }
    }

    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        BeamInfo *beaminfo    = params.beaminfo;
        Preprocess(params);

        if(beaminfo->SBPFlag == '*') {
            PRTFile
                << "\n______________________________\nUsing source beam pattern file\n";
        }
        // LP: Not normally echoed if there isn't a file, but may want to echo
        // values modified by the user.

        PRTFile << "Number of source beam pattern points " << beaminfo->NSBPPts << "\n";
        PRTFile << "\n Angle (degrees)  Power (dB)\n" << std::setprecision(3);
        for(int32_t i = 0; i < beaminfo->NSBPPts; ++i) {
            real dB = FL(20.0) * STD::log10(beaminfo->SrcBmPat[2 * i + 1]);
            PRTFile << beaminfo->SrcBmPat[2 * i] << " " << dB << "\n";
        }
    }

    virtual void Preprocess(bhcParams<O3D, R3D> &params) const
    {
        BeamInfo *beaminfo = params.beaminfo;
        if(beaminfo->SBPIndB) {
            beaminfo->SBPIndB = false;
            // convert dB to linear scale
            for(int32_t i = 0; i < beaminfo->NSBPPts; ++i)
                beaminfo->SrcBmPat[i * 2 + 1]
                    = STD::pow(FL(10.0), beaminfo->SrcBmPat[i * 2 + 1] / FL(20.0));
        }
    }

    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        trackdeallocate(params, params.beaminfo->SrcBmPat);
    }
};

}} // namespace bhc::module
