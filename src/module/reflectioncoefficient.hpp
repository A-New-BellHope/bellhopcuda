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
 * Optionally read in reflection coefficient for Top or Bottom boundary
 *
 * flag set to 'F' if refl. coef. is to be read from a File
 */
template<bool O3D, bool R3D, bool ISTOP> class ReflectionCoefficient {
public:
    ReflectionCoefficient() {}
    virtual ~ReflectionCoefficient() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        ReflectionInfoTopBot *refltb = GetReflTopBot(params);
        refltb->r                    = nullptr;
    }

    virtual void Default(bhcParams<O3D, R3D> &params) const
    {
        ReflectionInfoTopBot *refltb = GetReflTopBot(params);
        // LP: Original version allocates a size 1 array but never initializes
        // NPts. However, none of this is actually read in Reflect() unless
        // IsFile(), so it isn't actually wrong.
        refltb->NPts = 0;
    }

    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        ReflectionInfoTopBot *refltb = GetReflTopBot(params);
        if(!IsFile(params)) {
            Default(params);
            return;
        }
        LDIFile BRCFile(GetInternal(params), GetInternal(params)->FileRoot + s_extension);
        if(!BRCFile.Good()) {
            PRTFile << s_RC << "File = " << GetInternal(params)->FileRoot << s_extension
                    << "\n";
            EXTERR(
                "ReadReflectionCoefficient: Unable to open %s Reflection "
                "Coefficient file",
                s_TopBottom);
        }
        LIST(BRCFile);
        BRCFile.Read(refltb->NPts);
        trackallocate(params, "reflection coefficients", refltb->r, refltb->NPts);

        LIST(BRCFile);
        for(int32_t itheta = 0; itheta < refltb->NPts; ++itheta) {
            BRCFile.Read(refltb->r[itheta].theta);
            BRCFile.Read(refltb->r[itheta].r);
            BRCFile.Read(refltb->r[itheta].phi);
            refltb->r[itheta].phi *= DegRad; // convert to radians
        }
        refl->AnglesInDegrees = true;
    }

    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        if constexpr(!ISTOP) {
            // Optionally read in internal reflection coefficient data
            if(params.Bdry->Bot.hs.Opt[0] == 'P') {
                EXTERR("Internal reflections not supported by BELLHOP and therefore "
                       "not supported by " BHC_PROGRAMNAME);
            }
        }
    }

    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        if(!IsFile(params)) return;
        PRTFile << "_____________________________________________________________________"
                   "_____\n\n";
        PRTFile << "Using tabulated " << s_topbottom << " reflection coef.\n";
        PRTFile << "Number of points in " << s_topbottom
                << " reflection coefficient = " << refl->bot.NPts << "\n";
    }

    virtual void PostProcess(const bhcParams<O3D, R3D> &params) const
    {
        if(!IsFile(params)) return;

        ReflectionInfoTopBot *refltb = GetReflTopBot(params);
        if(refltb->AnglesInDegrees) {
            refltb->AnglesInDegrees = false;
            for(int32_t itheta = 0; itheta < refltb->NPts; ++itheta) {
                refltb->r[itheta].phi *= DegRad; // convert to radians
            }
        }
    }

    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        ReflectionInfoTopBot *refltb = GetReflTopBot(params);
        trackdeallocate(params, refltb->r);
    }

private:
    ReflectionInfoTopBot *GetReflTopBot(bhcParams<O3D, R3D> &params) const
    {
        if constexpr(ISTOP)
            return params.refl->top;
        else
            return params.refl->bot;
    }
    bool IsFile(bhcParams<O3D, R3D> &params) const
    {
        char opt;
        if constexpr(ISTOP)
            opt = params.Bdry->Top.hs.Opt[1];
        else
            opt = params.Bdry->Bot.hs.Opt[0];
        return opt == 'F';
    }
    constexpr const char *s_topbottom = ISTOP ? "top   " : "bottom";
    constexpr const char *s_TopBottom = ISTOP ? "Top" : "Bottom";
    constexpr const char *s_extension = ISTOP ? ".trc" : ".brc";
    constexpr const char *s_RC        = ISTOP ? "TRC" : "BRC";
};

template<bool O3D, bool R3D> using TRC = ReflectionCoefficient<O3D, R3D, true>;
template<bool O3D, bool R3D> using BRC = ReflectionCoefficient<O3D, R3D, false>;

}} // namespace bhc::module
