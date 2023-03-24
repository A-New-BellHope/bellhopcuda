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
#include "../common_setup.hpp"
#include "paramsmodule.hpp"

namespace bhc { namespace module {

template<bool O3D, bool R3D> class SSP {
public:
    SSP() {}
    virtual ~SSP() {}

    virtual void Init(bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp = params.ssp;

        ssp->cMat  = nullptr;
        ssp->czMat = nullptr;
        ssp->Seg.r = nullptr;
        ssp->Seg.x = nullptr;
        ssp->Seg.y = nullptr;
        ssp->Seg.z = nullptr;
    }

    virtual void SetupPre(bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp = params.ssp;

        ssp->NPts = 2;
        ssp->z[0] = RL(0.0);
        ssp->z[1] = RL(5000.0);

        ssp->Nr = ssp->Nx = ssp->Ny = 0;
        ssp->Nz                     = ssp->NPts;

        ssp->dirty     = true;
        ssp->rangeInKm = false;
    }

    virtual void Read(
        bhcParams<O3D, R3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const
    {
        SSPStructure *ssp = params.ssp;

        LIST(ENVFile);
        ENVFile.Read(NPts);
        ENVFile.Read(Sigma);
        ENVFile.Read(params.Bdry->Bot.hs.Depth);

        if(ssp->Type == 'A') return;

        ssp->NPts = 0;

        while(true) {
            if(ssp->NPts >= MaxSSP) {
                EXTERR("ReadSSP: Number of SSP points exceeds limit");
                return;
            }

            LIST_WARNLINE(ENVFile);
            ENVFile.Read(ssp->z[ssp->NPts]);
            ENVFile.Read(RecycledHS.alphaR);
            ENVFile.Read(RecycledHS.betaR);
            ENVFile.Read(RecycledHS.rho);
            ENVFile.Read(RecycledHS.alphaI);
            ENVFile.Read(RecycledHS.betaI);

            ssp->alphaR[ssp->NPts] = RecycledHS.alphaR;
            ssp->betaR[ssp->NPts]  = RecycledHS.betaR;
            ssp->rho[ssp->NPts]    = RecycledHS.rho;
            ssp->alphaI[ssp->NPts] = RecycledHS.alphaI;
            ssp->betaI[ssp->NPts]  = RecycledHS.betaI;

            ++ssp->NPts;

            // Did we read the last point?
            // LP: FLT_EPSILON is not a typo
            if(std::abs(ssp->z[ssp->NPts - 1] - params.Bdry->Bot.hs.Depth)
               < FL(100.0) * FLT_EPSILON) {
                break;
            }
        }

        ssp->Nz = ssp->NPts;

        if(ssp->Type == 'Q') {
            // read in extra SSP data for 2D
            ssp->rangeInKm = true;

            LDIFile SSPFile(GetInternal(params), GetInternal(params)->FileRoot + ".ssp");
            LIST(SSPFile);
            SSPFile.Read(ssp->Nr);

            trackallocate(params, "quad SSP", ssp->cMat, ssp->NPts * ssp->Nr);
            trackallocate(
                params, "quad SSP derivatives", ssp->czMat, (ssp->NPts - 1) * ssp->Nr);
            trackallocate(params, "quad SSP ranges", ssp->Seg.r, ssp->Nr);

            LIST(SSPFile);
            SSPFile.Read(ssp->Seg.r, ssp->Nr);

            for(int32_t iz2 = 0; iz2 < ssp->NPts; ++iz2) {
                LIST(SSPFile);
                SSPFile.Read(&ssp->cMat[iz2 * ssp->Nr], ssp->Nr);
            }
        } else if(ssp->Type == 'H') {
            // read in extra SSP data for 3D
            ssp->rangeInKm = true;

            // Read the 3D SSP matrix
            LDIFile SSPFile(GetInternal(params), GetInternal(params)->FileRoot + ".ssp");

            // x coordinates
            LIST(SSPFile);
            SSPFile.Read(ssp->Nx);
            trackallocate(params, "hexahedral SSP grid", ssp->Seg.x, ssp->Nx);
            LIST(SSPFile);
            SSPFile.Read(ssp->Seg.x, ssp->Nx);

            // y coordinates
            LIST(SSPFile);
            SSPFile.Read(ssp->Ny);
            trackallocate(params, "hexahedral SSP grid", ssp->Seg.y, ssp->Ny);
            LIST(SSPFile);
            SSPFile.Read(ssp->Seg.y, ssp->Ny);

            // z coordinates
            LIST(SSPFile);
            SSPFile.Read(ssp->Nz);
            trackallocate(params, "hexahedral SSP grid", ssp->Seg.z, ssp->Nz);
            LIST(SSPFile);
            SSPFile.Read(ssp->Seg.z, ssp->Nz);

            trackallocate(
                params, "hexahedral SSP values", ssp->cMat, ssp->Nx * ssp->Ny * ssp->Nz);
            trackallocate(
                params, "hexahedral SSP derivatives", ssp->czMat,
                ssp->Nx * ssp->Ny * (ssp->Nz - 1));

            for(int32_t iz2 = 0; iz2 < ssp->Nz; ++iz2) {
                for(int32_t iy2 = 0; iy2 < ssp->Ny; ++iy2) {
                    LIST(SSPFile);
                    for(int32_t ix2 = 0; ix2 < ssp->Nx; ++ix2) {
                        SSPFile.Read(ssp->cMat[((ix2)*ssp->Ny + iy2) * ssp->Nz + iz2]);
                    }
                }
            }
        }
    }

    virtual void SetupPost(bhcParams<O3D, R3D> &params) const override
    {
        // Depth of top boundary is taken from first SSP point
        // LP: Must be here, used later in env setup.
        Bdry->Top.hs.Depth = params.ssp->z[0];
        // [mbp:] bottom depth should perhaps be set the same way?
    }

    virtual void Validate(const bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp = params.ssp;

        if(ssp->NPts > MaxSSP) {
            EXTERR("ReadSSP: Number of SSP points exceeds limit");
        } else if(ssp->NPts < 2) {
            EXTERR("ReadSSP: The SSP must have at least 2 points");
        }

        if(ssp->NPts != ssp->Nz) {
            EXTERR("ssp->Npts / ssp->Nz have not been set up correctly");
        }

        if(ssp->z[0] != params.Bdry->Top.hs.Depth) {
            EXTERR("Ocean surface (Bdry->Top.hs.Depth) must be at first SSP depth "
                   "(ssp->z[0])");
        }
        if(std::abs(ssp->z[ssp->NPts - 1] - params.Bdry->Bot.hs.Depth)
           >= FL(100.0) * FLT_EPSILON) {
            EXTERR("Ocean bottom (Bdry->Bot.hs.Depth) must be at last SSP depth "
                   "(ssp->z[ssp->NPts-1])");
        }

        switch(ssp->Type) {
        case 'N': break;
        case 'C': break;
        case 'S': break;
        case 'P':
            if constexpr(O3D) {
#ifdef BHC_LIMIT_FEATURES
                EXTERR("PreprocessSSP: PCHIP is not supported in BELLHOP3D in "
                       "3D or Nx2D mode, but can be supported in " BHC_PROGRAMNAME
                       "if you turn off BHC_LIMIT_FEATURES");
#else
                EXTWARN("PreprocessSSP: warning: PCHIP not supported in BELLHOP3D in "
                        "3D or Nx2D mode, but supported in " BHC_PROGRAMNAME);
#endif
            }
        case 'Q':
            if constexpr(O3D) {
                EXTERR("PreprocessSSP: 2D profile not supported in 3D or Nx2D mode");
            }
            if(ssp->Nr < 2) {
                PRTFile << "ssp: Quad: You must have a least two profiles in your 2D SSP "
                           "field\n";
                std::abort();
            }
            break;
        case 'H':
            if constexpr(!O3D) {
                EXTERR("PreprocessSSP: 3D profile not supported in 2D mode");
            }
            if(ssp->Nx < 2 || ssp->Ny < 2 || ssp->Nz < 2) {
                EXTERR("ssp: Hexahedral: You must have at least two points in x, y, z "
                       "directions in your 3D SSP field");
            }
            if(ssp->Nz >= MaxSSP) {
                EXTERR("ssp: Hexahedral: Number of SSP points in Z exceeds limit");
            }
            break;
        case 'A': break;
        default: EXTERR("PreprocessSSP: Invalid profile option %c", ssp->Type);
        }

        if(!monotonic(ssp->NPts, ssp->z)) {
            EXTERR("PreprocessSSP: The depths in the SSP must be monotone increasing ");
        }
    }

    virtual void Echo(const bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp     = params.ssp;
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

        PRTFile << "\n  Depth = " << std::setw(10) << std::setprecision(2)
                << params.Bdry->Bot.hs.Depth << "  m\n";

        if(ssp->Type == 'A') {
            PRTFile << "Analytic SSP option\n";
            return;
        }

        PRTFile << "\nSound speed profile:\n";
        PRTFile << "      z         alphaR      betaR     rho        alphaI     betaI\n";
        PRTFile << "     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)\n";

        for(int32_t i = 0; i < ssp->NPts; ++i) {
            PRTFile << std::setprecision(2) << ssp->z[ssp->NPts] << " ";
            PRTFile << ssp->alphaR << " ";
            PRTFile << ssp->betaR << " ";
            PRTFile << ssp->rho << " ";
            PRTFile << std::setprecision(4) << ssp->alphaI << " ";
            PRTFile << ssp->betaI << "\n";
        }

        switch(ssp->Type) {
        case 'Q':
            PRTFile << "_________________________________________________________________"
                       "________"
                       "_\n\n";
            PRTFile << "Using range-dependent sound speed\n";

            PRTFile << "Number of SSP ranges = " << ssp->Nr << "\n";
            PRTFile << "\nProfile ranges (km):\n" << std::setprecision(2);
            for(int32_t i = 0; i < ssp->Nr; ++i) PRTFile << ssp->Seg.r[i] << " ";
            PRTFile << "\n";

            PRTFile << "\nSound speed matrix:\n";
            PRTFile << " Depth (m )     Soundspeed (m/s)\n";
            for(int32_t iz2 = 0; iz2 < ssp->NPts; ++iz2) {
                // PRTFile << "iSeg.z depth = " << std::setprecision(2) << ssp->z[iz2] <<
                // " m\n";
                PRTFile << std::setprecision(2) << ssp->z[iz2] << " ";
                for(int32_t i = 0; i < ssp->Nr; ++i)
                    PRTFile << ssp->cMat[iz2 * ssp->Nr + i] << " ";
                PRTFile << "\n";
            }
            break;
        case 'H':
            PRTFile << "\nReading sound speed profile from file\n";
            PRTFile << "\nNumber of points in x = " << ssp->Nx << "\n";
            PRTFile << "\nNumber of points in y = " << ssp->Ny << "\n";
            PRTFile << "\nNumber of points in z = " << ssp->Nz << "\n";
            PRTFile << "\n";
            break;
        }
    }

    virtual void Preprocess(bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp = params.ssp;

        if(ssp->rangeInKm) {
            // convert km to m
            if(ssp->Type == 'Q') {
                ToMeters2(ssp->Nr, ssp->Seg.r);
            } else if(ssp->Type == 'H') {
                ToMeters2(ssp->Nx, ssp->Seg.x);
                ToMeters2(ssp->Ny, ssp->Seg.y);
            }
            ssp->rangeInKm = false;
        }

        if(!ssp->dirty) return;
        ssp->dirty = false;

        if(ssp->Type == 'H') {
            // calculate cz
            for(int32_t iSegxt = 0; iSegxt < ssp->Nx; ++iSegxt) {
                for(int32_t iSegyt = 0; iSegyt < ssp->Ny; ++iSegyt) {
                    for(int32_t iz2 = 1; iz2 < ssp->Nz; ++iz2) {
                        ssp->czMat[((iSegxt)*ssp->Ny + iSegyt) * (ssp->Nz - 1) + iz2 - 1]
                            = (ssp->cMat[((iSegxt)*ssp->Ny + iSegyt) * ssp->Nz + iz2]
                               - ssp->cMat
                                     [((iSegxt)*ssp->Ny + iSegyt) * ssp->Nz + iz2 - 1])
                            / (ssp->Seg.z[iz2] - ssp->Seg.z[iz2 - 1]);
                    }
                }
            }

            // over-ride the SSP%z values read in from the environmental file with these
            // new values
            ssp->NPts = ssp->Nz;
            for(int32_t iz = 0; iz < ssp->Nz; ++iz) {
                ssp->z[iz] = ssp->Seg.z[iz];
                // LP: These are not well-defined, make sure they're not used
                ssp->c[iz]  = NAN;
                ssp->cz[iz] = NAN;
            }

            // LP: ssp->c and ssp->cz are not well-defined in hexahedral mode, and
            // if the number of depths is changed (ssp->Nz vs. ssp->NPts), computing
            // them may read uninitialized data.
            return;
        }

        for(int32_t iz = 0; iz < ssp->NPts; ++iz) {
            ssp->c[iz] = crci(
                params, ssp->z[iz], ssp->alphaR[iz], ssp->alphaI[iz], ssp->AttenUnit);

            // compute gradient, cz
            if(iz > 0)
                ssp->cz[iz - 1] = (ssp->c[iz] - ssp->c[iz - 1])
                    / (ssp->z[iz] - ssp->z[iz - 1]);
        }
        // LP: Gradient at last point is uninitialized.
        ssp->cz[ssp->NPts - 1] = cpx(NAN, NAN);

        switch(ssp->Type) {
        case 'N': // N2-linear profile option
            for(int32_t i = 0; i < ssp->NPts; ++i) ssp->n2[i] = FL(1.0) / SQ(ssp->c[i]);

            // compute gradient, n2z
            for(int32_t iz = 1; iz < ssp->NPts; ++iz) {
                ssp->n2z[iz - 1] = (ssp->n2[iz] - ssp->n2[iz - 1])
                    / (ssp->z[iz] - ssp->z[iz - 1]);
            }
            break;
        case 'S': { // Cubic spline profile option
            for(int32_t i = 0; i < ssp->NPts; ++i) ssp->cSpline[0][i] = ssp->c[i];

            // Compute spline coefs
            int32_t iBCBeg = 0;
            int32_t iBCEnd = 0;
            cSpline(
                ssp->z, ssp->cSpline[0], ssp->cSpline[1], ssp->cSpline[2],
                ssp->cSpline[3], ssp->NPts, iBCBeg, iBCEnd, ssp->NPts);
        } break;
        case 'P': // monotone PCHIP ACS profile option
            //                                                               2      3
            // compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
            //
            pchip(
                ssp->z, ssp->c, ssp->NPts, ssp->cCoef[0], ssp->cCoef[1], ssp->cCoef[2],
                ssp->cCoef[3], ssp->CSWork[0], ssp->CSWork[1], ssp->CSWork[2],
                ssp->CSWork[3]);
            break;
        case 'Q':
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
            break;
        }
    }

    virtual void Finalize(bhcParams<O3D, R3D> &params) const
    {
        SSPStructure *ssp = params.ssp;

        trackdeallocate(params, ssp->cMat);
        trackdeallocate(params, ssp->czMat);
        trackdeallocate(params, ssp->Seg.r);
        trackdeallocate(params, ssp->Seg.x);
        trackdeallocate(params, ssp->Seg.y);
        trackdeallocate(params, ssp->Seg.z);
    }
};

}} // namespace bhc::module
