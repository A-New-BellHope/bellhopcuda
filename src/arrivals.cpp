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
#include "arrivals.hpp"
#include "runtype.hpp"

namespace bhc {

template<bool O3D, bool R3D> void FinalizeArrivalsMode(
    const ArrInfo *arrinfo, const Position *Pos, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, std::string FileRoot)
{
    // LP: originally most of OpenOutputFiles
    bool isAscii;
    LDOFile AARRFile;
    UnformattedOFile BARRFile;
    switch(Beam->RunType[0]) {
    case 'A': // arrivals calculation, ascii
        isAscii = true;

        AARRFile.open(FileRoot + ".arr");
        AARRFile << (O3D ? "3D" : "2D") << '\n';
        AARRFile << freqinfo->freq0 << '\n';

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

        BARRFile.open(FileRoot + ".arr");
        BARRFile.rec();
        BARRFile.write((O3D ? "'3D'" : "'2D'"), 4);
        BARRFile.rec();
        BARRFile.write((float)freqinfo->freq0);

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
    default:
        GlobalLog("FinalizeArrivalsMode called while not in arrivals mode\n");
        bail();
        return;
    }
    // LP: originally most of WriteArrivals[ASCII/Binary][3D]
    for(int32_t isz = 0; isz < Pos->NSz; ++isz) {
        for(int32_t isx = 0; isx < Pos->NSx; ++isx) {
            for(int32_t isy = 0; isy < Pos->NSy; ++isy) {
                // LP: Maximum number of arrivals
                int32_t maxn = 0;
                for(int32_t itheta = 0; itheta < Pos->Ntheta; ++itheta) {
                    for(int32_t iz = 0; iz < Pos->NRz_per_range; ++iz) {
                        for(int32_t ir = 0; ir < Pos->NRr; ++ir) {
                            size_t base
                                = GetFieldAddr(isx, isy, isz, itheta, iz, ir, Pos);
                            if(arrinfo->NArr[base] > maxn) maxn = arrinfo->NArr[base];
                        }
                    }
                }
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
                            float factor;
                            if constexpr(R3D) {
                                factor = FL(1.0);
                            } else {
                                if(!O3D && IsLineSource(Beam)) {
                                    factor = FL(4.0) * STD::sqrt(REAL_PI);
                                } else if(Pos->Rr[ir] == FL(0.0)) {
                                    factor = FL(1e5); // avoid /0 at origin
                                } else {
                                    factor = FL(1.0)
                                        / STD::sqrt(Pos->Rr[ir]); // cyl. spreading
                                }
                            }

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
                                    AARRFile << factor * arr->a;
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
                                    BARRFile.write(factor * arr->a);
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
template void FinalizeArrivalsMode<false, false>(
    const ArrInfo *arrinfo, const Position *Pos, const FreqInfo *freqinfo,
    const BeamStructure<false> *Beam, std::string FileRoot);
#endif
#if BHC_ENABLE_NX2D
template void FinalizeArrivalsMode<true, false>(
    const ArrInfo *arrinfo, const Position *Pos, const FreqInfo *freqinfo,
    const BeamStructure<true> *Beam, std::string FileRoot);
#endif
#if BHC_ENABLE_3D
template void FinalizeArrivalsMode<true, true>(
    const ArrInfo *arrinfo, const Position *Pos, const FreqInfo *freqinfo,
    const BeamStructure<true> *Beam, std::string FileRoot);
#endif

} // namespace bhc
