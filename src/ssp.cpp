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
#include "ssp.hpp"
#include "curves.hpp"
#include "boundary.hpp"

namespace bhc {

void Initn2Linear(SSPStructure *ssp)
{
    for(int32_t i = 0; i < ssp->NPts; ++i) ssp->n2[i] = FL(1.0) / SQ(ssp->c[i]);

    // compute gradient, n2z
    for(int32_t iz = 1; iz < ssp->NPts; ++iz) {
        ssp->n2z[iz - 1] = (ssp->n2[iz] - ssp->n2[iz - 1])
            / (ssp->z[iz] - ssp->z[iz - 1]);
    }
}

void InitcPCHIP(SSPStructure *ssp)
{
    //                                                               2      3
    // compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
    //
    pchip(
        ssp->z, ssp->c, ssp->NPts, ssp->cCoef[0], ssp->cCoef[1], ssp->cCoef[2],
        ssp->cCoef[3], ssp->CSWork[0], ssp->CSWork[1], ssp->CSWork[2], ssp->CSWork[3]);
}

void InitcCubic(SSPStructure *ssp)
{
    for(int32_t i = 0; i < ssp->NPts; ++i) ssp->cSpline[0][i] = ssp->c[i];

    // Compute spline coefs
    int32_t iBCBeg = 0;
    int32_t iBCEnd = 0;
    cSpline(
        ssp->z, ssp->cSpline[0], ssp->cSpline[1], ssp->cSpline[2], ssp->cSpline[3],
        ssp->NPts, iBCBeg, iBCEnd, ssp->NPts);
}

template<bool O3D, bool R3D> void ReadQuad(bhcParams<O3D, R3D> &params)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    SSPStructure *ssp     = params.ssp;

    // Read the 2D SSP matrix
    PRTFile << "_________________________________________________________________________"
               "_\n\n";
    PRTFile << "Using range-dependent sound speed\n";

    LDIFile SSPFile(GetInternal(params), GetInternal(params)->FileRoot + ".ssp");
    LIST(SSPFile);
    SSPFile.Read(ssp->Nr);
    PRTFile << "Number of SSP ranges = " << ssp->Nr << "\n";

    if(ssp->Nr < 2) {
        PRTFile << "ssp: Quad: You must have a least two profiles in your 2D SSP field\n";
        std::abort();
    }

    trackallocate(params, "quad SSP", ssp->cMat, ssp->NPts * ssp->Nr);
    trackallocate(params, "quad SSP derivatives", ssp->czMat, (ssp->NPts - 1) * ssp->Nr);
    trackallocate(params, "quad SSP ranges", ssp->Seg.r, ssp->Nr);

    LIST(SSPFile);
    SSPFile.Read(ssp->Seg.r, ssp->Nr);
    PRTFile << "\nProfile ranges (km):\n" << std::setprecision(2);
    for(int32_t i = 0; i < ssp->Nr; ++i) PRTFile << ssp->Seg.r[i] << " ";
    PRTFile << "\n";

    for(int32_t i = 0; i < ssp->Nr; ++i) ssp->Seg.r[i] *= FL(1000.0); // convert km to m

    PRTFile << "\nSound speed matrix:\n";
    PRTFile << " Depth (m )     Soundspeed (m/s)\n";
    for(int32_t iz2 = 0; iz2 < ssp->NPts; ++iz2) {
        LIST(SSPFile);
        SSPFile.Read(&ssp->cMat[iz2 * ssp->Nr], ssp->Nr);
        // PRTFile << "iSeg.z depth = " << std::setprecision(2) << ssp->z[iz2] << " m\n";
        PRTFile << std::setprecision(2) << ssp->z[iz2] << " ";
        for(int32_t i = 0; i < ssp->Nr; ++i)
            PRTFile << ssp->cMat[iz2 * ssp->Nr + i] << " ";
        PRTFile << "\n";
    }
}

void InitQuad(SSPStructure *ssp)
{
    // calculate cz
    for(int32_t iSegt = 0; iSegt < ssp->Nr; ++iSegt) {
        for(int32_t iz2 = 1; iz2 < ssp->NPts; ++iz2) {
            real delta_z = ssp->z[iz2] - ssp->z[iz2 - 1];
            ssp->czMat[(iz2 - 1) * ssp->Nr + iSegt]
                = (ssp->cMat[iz2 * ssp->Nr + iSegt]
                   - ssp->cMat[(iz2 - 1) * ssp->Nr + iSegt])
                / delta_z;
        }
    }
}

template<bool O3D, bool R3D> void ReadHexahedral(bhcParams<O3D, R3D> &params)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
    SSPStructure *ssp     = params.ssp;

    // Read the 3D SSP matrix
    PRTFile << "\nReading sound speed profile from file\n";
    LDIFile SSPFile(GetInternal(params), GetInternal(params)->FileRoot + ".ssp");

    // x coordinates
    LIST(SSPFile);
    SSPFile.Read(ssp->Nx);
    PRTFile << "\nNumber of points in x = " << ssp->Nx << "\n";
    trackallocate(params, "hexahedral SSP grid", ssp->Seg.x, ssp->Nx);
    LIST(SSPFile);
    SSPFile.Read(ssp->Seg.x, ssp->Nx);

    // y coordinates
    LIST(SSPFile);
    SSPFile.Read(ssp->Ny);
    PRTFile << "\nNumber of points in y = " << ssp->Ny << "\n";
    trackallocate(params, "hexahedral SSP grid", ssp->Seg.y, ssp->Ny);
    LIST(SSPFile);
    SSPFile.Read(ssp->Seg.y, ssp->Ny);

    // z coordinates
    LIST(SSPFile);
    SSPFile.Read(ssp->Nz);
    PRTFile << "\nNumber of points in z = " << ssp->Nz << "\n";
    trackallocate(params, "hexahedral SSP grid", ssp->Seg.z, ssp->Nz);
    LIST(SSPFile);
    SSPFile.Read(ssp->Seg.z, ssp->Nz);

    if(ssp->Nx < 2 || ssp->Ny < 2 || ssp->Nz < 2) {
        EXTERR("ssp: Hexahedral: You must have at least two points in x, y, z "
               "directions in your 3D SSP field");
    }

    if(ssp->Nz >= MaxSSP) {
        EXTERR("ssp: Hexahedral: Number of SSP points in Z exceeds limit");
    }

    trackallocate(
        params, "hexahedral SSP values", ssp->cMat, ssp->Nx * ssp->Ny * ssp->Nz);
    trackallocate(
        params, "hexahedral SSP derivatives", ssp->czMat,
        ssp->Nx * ssp->Ny * (ssp->Nz - 1));

    PRTFile << "\n";
    for(int32_t iz2 = 0; iz2 < ssp->Nz; ++iz2) {
        for(int32_t iy2 = 0; iy2 < ssp->Ny; ++iy2) {
            LIST(SSPFile);
            for(int32_t ix2 = 0; ix2 < ssp->Nx; ++ix2) {
                SSPFile.Read(ssp->cMat[((ix2)*ssp->Ny + iy2) * ssp->Nz + iz2]);
            }
        }
    }

    // convert km to m
    for(int32_t ix1 = 0; ix1 < ssp->Nx; ++ix1) ssp->Seg.x[ix1] *= FL(1000.0);
    for(int32_t iy1 = 0; iy1 < ssp->Ny; ++iy1) ssp->Seg.y[iy1] *= FL(1000.0);
}

void InitHexahedral(SSPStructure *ssp)
{
    // calculate cz
    for(int32_t iSegxt = 0; iSegxt < ssp->Nx; ++iSegxt) {
        for(int32_t iSegyt = 0; iSegyt < ssp->Ny; ++iSegyt) {
            for(int32_t iz2 = 1; iz2 < ssp->Nz; ++iz2) {
                ssp->czMat[((iSegxt)*ssp->Ny + iSegyt) * (ssp->Nz - 1) + iz2 - 1]
                    = (ssp->cMat[((iSegxt)*ssp->Ny + iSegyt) * ssp->Nz + iz2]
                       - ssp->cMat[((iSegxt)*ssp->Ny + iSegyt) * ssp->Nz + iz2 - 1])
                    / (ssp->Seg.z[iz2] - ssp->Seg.z[iz2 - 1]);
            }
        }
    }

    // over-ride the SSP%z values read in from the environmental file with these new
    // values
    ssp->NPts = ssp->Nz;
    for(int32_t iz = 0; iz < ssp->Nz; ++iz) {
        ssp->z[iz] = ssp->Seg.z[iz];
        // LP: These are not well-defined, make sure they're not used
        ssp->c[iz]  = DEBUG_LARGEVAL;
        ssp->cz[iz] = DEBUG_LARGEVAL;
    }
}

/**
 * Update the SSP parameters. Safe to call multiple times with flags.
 * Be sure to flag ssp->dirty if you change the SSP externally.
 */
template<bool O3D, bool R3D> void UpdateSSP(bhcParams<O3D, R3D> &params)
{
    SSPStructure *ssp = params.ssp;

    if(!ssp->dirty) return;
    ssp->dirty = false;

    if(ssp->Type == 'H') {
        if constexpr(!O3D) { EXTERR("UpdateSSP: 3D profile not supported in 2D mode"); }
        // LP: ssp->c and ssp->cz are not well-defined in hexahedral mode, and
        // if the number of depths is changed (ssp->Nz vs. ssp->NPts), computing
        // them may read uninitialized data.
        InitHexahedral(ssp);
    } else {
        for(int32_t iz = 0; iz < ssp->NPts; ++iz) {
            ssp->c[iz] = crci(
                params, ssp->z[iz], ssp->alphaR[iz], ssp->alphaI[iz], ssp->AttenUnit);

            // verify that the depths are monotone increasing
            if(iz > 0) {
                if(ssp->z[iz] <= ssp->z[iz - 1]) {
                    EXTERR(
                        "UpdateSSP: The depths in the SSP must be monotone increasing "
                        "(%f)",
                        ssp->z[iz]);
                }
            }

            // compute gradient, cz
            if(iz > 0)
                ssp->cz[iz - 1] = (ssp->c[iz] - ssp->c[iz - 1])
                    / (ssp->z[iz] - ssp->z[iz - 1]);
        }
        // LP: Gradient at last point is uninitialized.
        ssp->cz[ssp->NPts - 1] = cpx(FL(5.5555555e30), FL(-3.3333333e29)); // LP:
                                                                           // debugging

        switch(ssp->Type) {
        case 'N': // N2-linear profile option
            Initn2Linear(ssp);
            break;
        case 'C': // C-linear profile option
            // nothing to do
            break;
        case 'S': // Cubic spline profile option
            InitcCubic(ssp);
            break;
        case 'P': // monotone PCHIP ACS profile option
            if constexpr(O3D) {
#ifdef BHC_LIMIT_FEATURES
                EXTERR("UpdateSSP: PCHIP is not supported in BELLHOP3D in "
                       "3D or Nx2D mode, but can be supported in " BHC_PROGRAMNAME
                       "if you turn off BHC_LIMIT_FEATURES");
#else
                EXTWARN("UpdateSSP: warning: PCHIP not supported in BELLHOP3D in "
                        "3D or Nx2D mode, but supported in " BHC_PROGRAMNAME);
#endif
            }
            InitcPCHIP(ssp);
            break;
        case 'Q':
            if constexpr(O3D) {
                EXTERR("UpdateSSP: 2D profile not supported in 3D or Nx2D mode");
            }
            InitQuad(ssp);
            break;
        default: EXTERR("UpdateSSP: Invalid profile option %c", ssp->Type);
        }
    }
}

#if BHC_ENABLE_2D
template void UpdateSSP<false, false>(bhcParams<false, false> &params);
#endif
#if BHC_ENABLE_NX2D
template void UpdateSSP<true, false>(bhcParams<true, false> &params);
#endif
#if BHC_ENABLE_3D
template void UpdateSSP<true, true>(bhcParams<true, true> &params);
#endif

/**
 * reads the SSP data from the environmental file and convert to Nepers/m
 */
template<bool O3D, bool R3D> inline void ReadSSP(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    PRTFile << "\nSound speed profile:\n";
    PRTFile << "      z         alphaR      betaR     rho        alphaI     betaI\n";
    PRTFile << "     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)\n";

    params.ssp->NPts = 1;

    for(int32_t iz = 0; iz < MaxSSP; ++iz) {
        LIST_WARNLINE(ENVFile);
        ENVFile.Read(params.ssp->z[iz]);
        ENVFile.Read(RecycledHS.alphaR);
        ENVFile.Read(RecycledHS.betaR);
        ENVFile.Read(RecycledHS.rho);
        ENVFile.Read(RecycledHS.alphaI);
        ENVFile.Read(RecycledHS.betaI);

        PRTFile << std::setprecision(2) << params.ssp->z[iz] << " " << RecycledHS.alphaR
                << " " << RecycledHS.betaR << " " << RecycledHS.rho << " "
                << std::setprecision(4) << RecycledHS.alphaI << " " << RecycledHS.betaI
                << "\n";
        params.ssp->rho[iz]    = RecycledHS.rho;
        params.ssp->alphaR[iz] = RecycledHS.alphaR;
        params.ssp->alphaI[iz] = RecycledHS.alphaI;

        // Did we read the last point?
        // LP: FLT_EPSILON is not a typo
        if(std::abs(params.ssp->z[iz] - params.Bdry->Bot.hs.Depth)
           < FL(100.0) * FLT_EPSILON) {
            params.ssp->Nz = params.ssp->NPts;
            if(params.ssp->NPts == 1) {
                EXTERR("ReadSSP: The SSP must have at least 2 points");
            }

            if(params.ssp->Type == 'Q') {
                // read in extra SSP data for 2D
                ReadQuad(params);
            } else if(params.ssp->Type == 'H') {
                // read in extra SSP data for 3D
                ReadHexahedral(params);
            }

            return;
        }

        ++params.ssp->NPts;
    }

    // Fall through means too many points in the profile
    EXTERR("ReadSSP: Number of SSP points exceeds limit");
}

template<bool O3D, bool R3D> void InitializeSSP(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS)
{
    if(params.ssp->Type == 'A') {
        // nothing to do for analytic
        return;
    }

    ReadSSP(params, ENVFile, RecycledHS);
    params.ssp->dirty = true;
}

#if BHC_ENABLE_2D
template void InitializeSSP<false, false>(
    bhcParams<false, false> &params, LDIFile &ENVFile, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_NX2D
template void InitializeSSP<true, false>(
    bhcParams<true, false> &params, LDIFile &ENVFile, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_3D
template void InitializeSSP<true, true>(
    bhcParams<true, true> &params, LDIFile &ENVFile, HSInfo &RecycledHS);
#endif

} // namespace bhc
