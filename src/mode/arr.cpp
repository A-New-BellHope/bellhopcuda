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
#include "arr.hpp"
#include "../common_run.hpp"

namespace bhc { namespace mode {

template<bool O3D, bool R3D> void PostProcessArrivals(
    const bhcParams<O3D, R3D> &params, ArrInfo *arrinfo)
{
    const Position *Pos = params.Pos;
    for(int32_t isz = 0; isz < Pos->NSz; ++isz) {
        for(int32_t isx = 0; isx < Pos->NSx; ++isx) {
            for(int32_t isy = 0; isy < Pos->NSy; ++isy) {
                int32_t maxn = 0; // LP: Maximum number of arrivals for this source
                for(int32_t itheta = 0; itheta < Pos->Ntheta; ++itheta) {
                    for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
                        for(int32_t ir = 0; ir < Pos->NRr; ++ir) {
                            size_t base
                                = GetFieldAddr(isx, isy, isz, itheta, iz, ir, Pos);

                            int32_t narr = arrinfo->NArr[base];
                            if(narr > arrinfo->MaxNArr) {
                                // For multithreading / AllowMerging == false where this
                                // holds the total number of attempted arrivals,
                                // including those not written due to limited memory
                                arrinfo->NArr[base] = narr = arrinfo->MaxNArr;
                            }
                            maxn = bhc::max(maxn, narr);

                            float factor;
                            if constexpr(R3D) {
                                factor = FL(1.0);
                            } else {
                                if(!O3D && IsLineSource(params.Beam)) {
                                    factor = FL(4.0) * STD::sqrt(REAL_PI);
                                } else if(Pos->Rr[ir] == FL(0.0)) {
                                    // avoid /0 at origin
                                    factor = FL(1e5);
                                } else {
                                    // cyl. spreading
                                    factor = FL(1.0) / STD::sqrt(Pos->Rr[ir]);
                                }
                            }
                            for(int32_t iArr = 0; iArr < narr; ++iArr) {
                                arrinfo->Arr[base * arrinfo->MaxNArr + iArr].a *= factor;
                            }
                        }
                    }
                }
                arrinfo->MaxNPerSource[(isz * Pos->NSx + isx) * Pos->NSy + isy] = maxn;
            }
        }
    }
}

#if BHC_ENABLE_2D
template void PostProcessArrivals<false, false>(
    const bhcParams<false, false> &params, ArrInfo *arrinfo);
#endif
#if BHC_ENABLE_NX2D
template void PostProcessArrivals<true, false>(
    const bhcParams<true, false> &params, ArrInfo *arrinfo);
#endif
#if BHC_ENABLE_3D
template void PostProcessArrivals<true, true>(
    const bhcParams<true, true> &params, ArrInfo *arrinfo);
#endif

template<bool O3D, bool R3D> void WriteOutArrivals(
    const bhcParams<O3D, R3D> &params, const ArrInfo *arrinfo)
{
    const Position *Pos = params.Pos;

    // LP: originally most of OpenOutputFiles
    bool isAscii;
    LDOFile AARRFile;
    UnformattedOFile BARRFile(GetInternal(params));
    switch(params.Beam->RunType[0]) {
    case 'A': // arrivals calculation, ascii
        isAscii = true;

        AARRFile.open(GetInternal(params)->FileRoot + ".arr");
        AARRFile << (O3D ? "3D" : "2D") << '\n';
        AARRFile << params.freqinfo->freq0 << '\n';

        // write source locations
        if constexpr(O3D) {
            AARRFile << Pos->NSx;
            AARRFile.write(Pos->Sx, Pos->NSx);
            AARRFile << '\n';
            AARRFile << Pos->NSy;
            AARRFile.write(Pos->Sy, Pos->NSy);
            AARRFile << '\n';
        }
        AARRFile << Pos->NSz;
        AARRFile.write(Pos->Sz, Pos->NSz);
        AARRFile << '\n';

        // write receiver locations
        AARRFile << Pos->NRz;
        AARRFile.write(Pos->Rz, Pos->NRz);
        AARRFile << '\n';
        AARRFile << Pos->NRr;
        AARRFile.write(Pos->Rr, Pos->NRr);
        AARRFile << '\n';
        if constexpr(O3D) {
            AARRFile << Pos->Ntheta;
            AARRFile.write(Pos->theta, Pos->Ntheta);
            AARRFile << '\n';
        }
        break;
    case 'a': // arrivals calculation, binary
        isAscii = false;

        BARRFile.open(GetInternal(params)->FileRoot + ".arr");
        BARRFile.rec();
        BARRFile.write((O3D ? "'3D'" : "'2D'"), 4);
        BARRFile.rec();
        BARRFile.write((float)params.freqinfo->freq0);

        // write source locations
        if constexpr(O3D) {
            BARRFile.rec();
            BARRFile.write(Pos->NSx);
            BARRFile.write(Pos->Sx, Pos->NSx);
            BARRFile.rec();
            BARRFile.write(Pos->NSy);
            BARRFile.write(Pos->Sy, Pos->NSy);
        }
        BARRFile.rec();
        BARRFile.write(Pos->NSz);
        BARRFile.write(Pos->Sz, Pos->NSz);

        // write receiver locations
        BARRFile.rec();
        BARRFile.write(Pos->NRz);
        BARRFile.write(Pos->Rz, Pos->NRz);
        BARRFile.rec();
        BARRFile.write(Pos->NRr);
        BARRFile.write(Pos->Rr, Pos->NRr);
        if constexpr(O3D) {
            BARRFile.rec();
            BARRFile.write(Pos->Ntheta);
            BARRFile.write(Pos->theta, Pos->Ntheta);
        }
        break;
    default: EXTERR("FinalizeArrivalsMode called while not in arrivals mode");
    }
    // LP: originally most of WriteArrivals[ASCII/Binary][3D]
    for(int32_t isz = 0; isz < Pos->NSz; ++isz) {
        for(int32_t isx = 0; isx < Pos->NSx; ++isx) {
            for(int32_t isy = 0; isy < Pos->NSy; ++isy) {
                // LP: Maximum number of arrivals for this source
                int32_t maxn
                    = arrinfo->MaxNPerSource[(isz * Pos->NSx + isx) * Pos->NSy + isy];
                if(isAscii) {
                    AARRFile << maxn << '\n';
                } else {
                    BARRFile.rec();
                    BARRFile.write(maxn);
                }

                for(int32_t itheta = 0; itheta < Pos->Ntheta; ++itheta) {
                    for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
                        for(int32_t ir = 0; ir < Pos->NRr; ++ir) {
                            size_t base
                                = GetFieldAddr(isx, isy, isz, itheta, iz, ir, Pos);
                            int32_t narr = arrinfo->NArr[base];
                            if(isAscii) {
                                AARRFile << narr << '\n';
                            } else {
                                BARRFile.rec();
                                BARRFile.write(narr);
                            }

                            for(int32_t iArr = 0; iArr < narr; ++iArr) {
                                Arrival *arr
                                    = &arrinfo->Arr[base * arrinfo->MaxNArr + iArr];
                                // LP: Unnecessary inconsistent casting to float; see
                                // Fortran version readme.
                                if(isAscii) {
                                    // You can compress the output file a lot by putting
                                    // in an explicit format statement here ... However,
                                    // you'll need to make sure you keep adequate
                                    // precision
                                    AARRFile << arr->a;
                                    if constexpr(O3D) {
                                        AARRFile << RadDeg * arr->Phase;
                                    } else {
                                        AARRFile << (float)RadDeg * arr->Phase;
                                    }
                                    AARRFile << arr->delay.real() << arr->delay.imag()
                                             << arr->SrcDeclAngle;
                                    if constexpr(O3D) AARRFile << arr->SrcAzimAngle;
                                    AARRFile << arr->RcvrDeclAngle;
                                    if constexpr(O3D) AARRFile << arr->RcvrAzimAngle;
                                    AARRFile << arr->NTopBnc << arr->NBotBnc << '\n';
                                } else {
                                    BARRFile.rec();
                                    BARRFile.write(arr->a);
                                    BARRFile.write((float)(RadDeg * arr->Phase));
                                    BARRFile.write(arr->delay);
                                    BARRFile.write(arr->SrcDeclAngle);
                                    if constexpr(O3D) BARRFile.write(arr->SrcAzimAngle);
                                    BARRFile.write(arr->RcvrDeclAngle);
                                    if constexpr(O3D) BARRFile.write(arr->RcvrAzimAngle);
                                    BARRFile.write((float)arr->NTopBnc);
                                    BARRFile.write((float)arr->NBotBnc);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#if BHC_ENABLE_2D
template void WriteOutArrivals<false, false>(
    const bhcParams<false, false> &params, const ArrInfo *arrinfo);
#endif
#if BHC_ENABLE_NX2D
template void WriteOutArrivals<true, false>(
    const bhcParams<true, false> &params, const ArrInfo *arrinfo);
#endif
#if BHC_ENABLE_3D
template void WriteOutArrivals<true, true>(
    const bhcParams<true, true> &params, const ArrInfo *arrinfo);
#endif

}} // namespace bhc::mode
