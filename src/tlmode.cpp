/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
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
#include "tlmode.hpp"

namespace bhc {

/**
 * Write header to disk file
 * LP: of a SHDFile ("binary `shade' file (SHDFIL) that contains calculated
 * pressure fields")
 *
 * FileName: Name of the file (could be a shade file or a Green's function file)
 * Title: Arbitrary title
 * freq0: Nominal frequency [LP: now in freqinfo]
 * atten: stabilizing attenuation (for wavenumber integration only)
 * PlotType: If "TL", writes only first and last Sx and Sy [LP: never set to
 * "TL" in BELLHOP]
 */
void WriteHeader(
    DirectOFile &SHDFile, const std::string &FileName, const char (&Title)[80],
    float atten, const std::string &PlotType, const Position *Pos,
    const FreqInfo *freqinfo)
{
    bool isTL = (PlotType[0] == 'T' && PlotType[1] == 'L');

    int32_t LRecl = 84; // 4 for LRecl, 80 for Title
    LRecl = bhc::max(LRecl, 2 * freqinfo->Nfreq * (int32_t)sizeof(freqinfo->freqVec[0]));
    LRecl = bhc::max(LRecl, Pos->Ntheta * (int32_t)sizeof(Pos->theta[0]));
    if(!isTL) {
        LRecl = bhc::max(LRecl, Pos->NSx * (int32_t)sizeof(Pos->Sx[0]));
        LRecl = bhc::max(LRecl, Pos->NSy * (int32_t)sizeof(Pos->Sy[0]));
    }
    LRecl = bhc::max(LRecl, Pos->NSz * (int32_t)sizeof(Pos->Sz[0]));
    LRecl = bhc::max(LRecl, Pos->NRz * (int32_t)sizeof(Pos->Rz[0]));
    LRecl = bhc::max(LRecl, Pos->NRr * (int32_t)sizeof(cpxf));

    SHDFile.open(FileName, LRecl);
    if(!SHDFile.good()) {
        GlobalLog("Could not open SHDFile: %s\n", FileName.c_str());
        std::abort();
    }
    LRecl /= 4;
    SHDFile.rec(0);
    DOFWRITE(SHDFile, &LRecl, 4);
    DOFWRITE(SHDFile, std::string(Title), 80);
    SHDFile.rec(1);
    DOFWRITE(SHDFile, PlotType, 10);
    SHDFile.rec(2);
    DOFWRITEV(SHDFile, freqinfo->Nfreq);
    DOFWRITEV(SHDFile, Pos->Ntheta);
    DOFWRITEV(SHDFile, Pos->NSx);
    DOFWRITEV(SHDFile, Pos->NSy);
    DOFWRITEV(SHDFile, Pos->NSz);
    DOFWRITEV(SHDFile, Pos->NRz);
    DOFWRITEV(SHDFile, Pos->NRr);
    DOFWRITEV(SHDFile, (float)freqinfo->freq0);
    DOFWRITEV(SHDFile, atten);
    SHDFile.rec(3);
    DOFWRITE(SHDFile, freqinfo->freqVec, freqinfo->Nfreq * sizeof(freqinfo->freqVec[0]));
    SHDFile.rec(4);
    DOFWRITE(SHDFile, Pos->theta, Pos->Ntheta * sizeof(Pos->theta[0]));

    if(!isTL) {
        SHDFile.rec(5);
        DOFWRITE(SHDFile, Pos->Sx, Pos->NSx * sizeof(Pos->Sx[0]));
        SHDFile.rec(6);
        DOFWRITE(SHDFile, Pos->Sy, Pos->NSy * sizeof(Pos->Sy[0]));
    } else {
        SHDFile.rec(5);
        DOFWRITEV(SHDFile, Pos->Sx[0]);
        DOFWRITEV(SHDFile, Pos->Sx[Pos->NSx - 1]);
        SHDFile.rec(6);
        DOFWRITEV(SHDFile, Pos->Sy[0]);
        DOFWRITEV(SHDFile, Pos->Sy[Pos->NSy - 1]);
    }
    SHDFile.rec(7);
    DOFWRITE(SHDFile, Pos->Sz, Pos->NSz * sizeof(Pos->Sz[0]));

    SHDFile.rec(8);
    DOFWRITE(SHDFile, Pos->Rz, Pos->NRz * sizeof(Pos->Rz[0]));
    SHDFile.rec(9);
    DOFWRITE(SHDFile, Pos->Rr, Pos->NRr * sizeof(Pos->Rr[0]));
}

/**
 * LP: Write TL results
 */
template<bool O3D, bool R3D> void WriteOutTL(
    const bhcParams<O3D, R3D> &params, bhcOutputs<O3D, R3D> &outputs)
{
    real atten = FL(0.0);
    std::string PlotType;
    DirectOFile SHDFile;

    // following to set PlotType has already been done in READIN if that was used for
    // input
    PlotType = IsIrregularGrid(params.Beam) ? "irregular " : "rectilin  ";
    WriteHeader(
        SHDFile, GetInternal(params)->FileRoot + ".shd", params.Title, atten, PlotType,
        params.Pos, params.freqinfo);

    for(int32_t isz = 0; isz < params.Pos->NSz; ++isz) {
        for(int32_t isx = 0; isx < params.Pos->NSx; ++isx) {
            for(int32_t isy = 0; isy < params.Pos->NSy; ++isy) {
                SSPSegState iSeg;
                iSeg.r = iSeg.x = iSeg.y = iSeg.z = 0;
                VEC23<O3D> xs, tinit;
                SSPOutputs<O3D> o;
                char st = params.ssp->Type;
                if(st == 'N') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'N'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'C') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'C'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'S') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'S'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'P') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'P'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'Q') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'Q'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'H') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'H'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else if(st == 'A') {
                    o = RayStartNominalSSP<CfgSel<'C', 'G', 'A'>, O3D>(
                        isx, isy, isz, FL(0.0), iSeg, params.Pos, params.ssp, xs, tinit);
                } else {
                    GlobalLog("Invalid ssp->Type %c!", st);
                    bail();
                }
                cpx epsilon1, epsilon2;
                if constexpr(R3D) {
                    epsilon1 = PickEpsilon<O3D, R3D>(
                        FL(2.0) * REAL_PI * params.freqinfo->freq0, o.ccpx.real(),
                        o.gradc, FL(0.0), params.Angles->alpha.d, params.Beam);
                    epsilon2 = PickEpsilon<O3D, R3D>(
                        FL(2.0) * REAL_PI * params.freqinfo->freq0, o.ccpx.real(),
                        o.gradc, FL(0.0), params.Angles->beta.d, params.Beam);
                } else {
                    epsilon1 = epsilon2 = RL(0.0);
                }
                ScalePressure<O3D, R3D>(
                    params.Angles->alpha.d, params.Angles->beta.d, o.ccpx.real(),
                    epsilon1, epsilon2, params.Pos->Rr,
                    &outputs
                         .uAllSources[GetFieldAddr(isx, isy, isz, 0, 0, 0, params.Pos)],
                    params.Pos->Ntheta, params.Pos->NRz_per_range, params.Pos->NRr,
                    params.freqinfo->freq0, params.Beam);
                // LP: depth
                for(int32_t Irz1 = 0; Irz1 < params.Pos->NRz_per_range; ++Irz1) {
                    for(int32_t itheta = 0; itheta < params.Pos->Ntheta; ++itheta) {
                        // LP: This is in a different order from the field.
                        // Field: (largest) Z, X, Y, theta, depth, radius (smallest)
                        // File:  (largest) X, Y, theta, Z, depth, radius (smallest)
                        // (XYZ are source; theta depth radius are receiver)
                        // Also, for 3D, mbp iterates over theta inside depth, which does
                        // not match the order of writing to the file.
                        // clang-format off
                        size_t IRec = 10                     + ((((size_t)isx
                            * (size_t)params.Pos->NSy           + (size_t)isy)
                            * (size_t)params.Pos->Ntheta        + (size_t)itheta)
                            * (size_t)params.Pos->NSz           + (size_t)isz)
                            * (size_t)params.Pos->NRz_per_range + (size_t)Irz1;
                        // clang-format on
                        SHDFile.rec(IRec);
                        for(int32_t r = 0; r < params.Pos->NRr; ++r) {
                            cpxf v = outputs.uAllSources[GetFieldAddr(
                                isx, isy, isz, itheta, Irz1, r, params.Pos)];
                            DOFWRITEV(SHDFile, v);
                        }
                    }
                }
            }
        }
    }
}

#if BHC_ENABLE_2D
template void WriteOutTL<false, false>(
    const bhcParams<false, false> &params, bhcOutputs<false, false> &outputs);
#endif
#if BHC_ENABLE_NX2D
template void WriteOutTL<true, false>(
    const bhcParams<true, false> &params, bhcOutputs<true, false> &outputs);
#endif
#if BHC_ENABLE_3D
template void WriteOutTL<true, true>(
    const bhcParams<true, true> &params, bhcOutputs<true, true> &outputs);
#endif

} // namespace bhc
