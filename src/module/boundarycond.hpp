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
#include "paramsmodule.hpp"

namespace bhc { namespace module {

/**
 * LP: Formerly TopBot
 * Handles top and bottom boundary conditions
 * params.freqinfo->freq0: center / nominal frequency (wideband not supported)
 */
template<bool O3D, bool ISTOP> class BoundaryCond : public ParamsModule<O3D> {
public:
    BoundaryCond() {}
    virtual ~BoundaryCond() {}

    virtual void SetupPre(bhcParams<O3D> &params) const override
    {
        HSInfo &hs = GetBdry(params).hs;
        hs.cP = hs.cS = hs.rho = FL(0.0);
        params.fT              = RL(1.0e20);
    }

    virtual void Default(bhcParams<O3D> &) const override {}

    virtual void Read(
        bhcParams<O3D> &params, LDIFile &ENVFile, HSInfo &RecycledHS) const override
    {
        HSInfo &hs   = GetBdry(params).hs;
        HSExtra &hsx = GetBdry(params).hsx;

        // ****** Read in BC parameters depending on particular choice ******

        if(hs.bc == 'A') { // *** Half-space properties ***
            hsx.zTemp = FL(0.0);
            LIST(ENVFile);
            ENVFile.Read(hsx.zTemp);
            ENVFile.Read(RecycledHS.alphaR);
            ENVFile.Read(RecycledHS.betaR);
            ENVFile.Read(RecycledHS.rho);
            ENVFile.Read(RecycledHS.alphaI);
            ENVFile.Read(RecycledHS.betaI);
            hs.alphaR = RecycledHS.alphaR;
            hs.betaR  = RecycledHS.betaR;
            hs.rho    = RecycledHS.rho;
            hs.alphaI = RecycledHS.alphaI;
            hs.betaI  = RecycledHS.betaI;

            // dummy parameters for a layer with a general power law for attenuation
            // these are not in play because the AttenUnit for this is not allowed yet
            // freq0         = freq;
            params.fT = FL(1000.0);
        } else if(hs.bc == 'G') { // *** Grain size (formulas from UW-APL HF Handbook)
            LIST(ENVFile);
            ENVFile.Read(hsx.zTemp);
            ENVFile.Read(hsx.Mz);
            Preprocess(params);
            RecycledHS.alphaR = hs.alphaR;
            RecycledHS.alphaI = hs.alphaI;
            RecycledHS.rho    = hs.rho;
        }
    }

    virtual void Write(bhcParams<O3D> &params, LDOFile &ENVFile) const
    {
        HSInfo &hs   = GetBdry(params).hs;
        HSExtra &hsx = GetBdry(params).hsx;

        if(hs.bc == 'A') { // *** Half-space properties ***
            ENVFile << hsx.zTemp;
            ENVFile << hs.alphaR << hs.betaR << hs.rho << hs.alphaI << hs.betaI;
            ENVFile.write("! zTemp alphaR betaR rho alphaI betaI\n");
        } else if(hs.bc == 'G') { // *** Grain size (formulas from UW-APL HF Handbook)
            ENVFile << hsx.zTemp << hsx.Mz;
            ENVFile.write("! zTemp Mz\n");
        }
    }

    virtual void Validate(bhcParams<O3D> &params) const override
    {
        switch(GetBdry(params).hs.bc) {
        case 'V': break;
        case 'R': break;
        case 'A': break;
        case 'G': break;
        case 'F': break;
        case 'W': break;
        case 'P': break;
        default: EXTERR("TopBot: Unknown boundary condition type");
        }
    }

    static void WriteBCTag(char opt, LDOFile &ENVFile)
    {
        switch(opt) {
        case 'V': ENVFile.write("vacuum"); break;
        case 'R': ENVFile.write("rigid"); break;
        case 'A': ENVFile.write("acousto-elastic"); break;
        case 'G': ENVFile.write("grain size"); break;
        case 'F': ENVFile.write("file"); break;
        case 'W': ENVFile.write("write IFL"); break;
        case 'P': ENVFile.write("precalculated IFL"); break;
        default: ENVFile.write("error!"); break;
        }
    }

    virtual void Echo(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;
        HSInfo &hs            = GetBdry(params).hs;
        HSExtra &hsx          = GetBdry(params).hsx;

        Preprocess(params);

        // Echo to PRTFile user's choice of boundary condition
        switch(hs.bc) {
        case 'V': PRTFile << "    VACUUM\n"; break;
        case 'R': PRTFile << "    Perfectly RIGID\n"; break;
        case 'A': PRTFile << "    ACOUSTO-ELASTIC half-space\n"; break;
        case 'G': PRTFile << "    Grain size to define half-space\n"; break;
        case 'F': PRTFile << "    FILE used for reflection loss\n"; break;
        case 'W': PRTFile << "    Writing an IFL file\n"; break;
        case 'P': PRTFile << "    reading PRECALCULATED IFL\n"; break;
        }

        if constexpr(ISTOP) {
            if(hs.bc == 'A') {
                PRTFile << "      z         alphaR      betaR     rho        alphaI     "
                           "betaI\n";
                PRTFile << "     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     "
                           "(m/s)\n";
            }
        }

        if(hs.bc == 'A') { // *** Half-space properties ***
            PRTFile << std::setprecision(2) << std::setw(10) << hsx.zTemp << " ";
            PRTFile << std::setw(10) << hs.alphaR << " ";
            PRTFile << std::setw(10) << hs.betaR << " ";
            PRTFile << std::setw(6) << hs.rho << " ";
            PRTFile << std::setprecision(4) << std::setw(10) << hs.alphaI << " ";
            PRTFile << std::setw(10) << hs.betaI << "\n";
        } else if(hs.bc == 'G') { // *** Grain size (formulas from UW-APL HF Handbook)
            PRTFile << std::setprecision(2) << std::setw(10) << hsx.zTemp << " ";
            PRTFile << std::setw(10) << hsx.Mz << "\n";
        }
    }

    virtual void Preprocess(bhcParams<O3D> &params) const override
    {
        HSInfo &hs   = GetBdry(params).hs;
        HSExtra &hsx = GetBdry(params).hsx;

        if(hs.bc == 'A') { // *** Half-space properties ***
            hs.cP = crci(params, hsx.zTemp, hs.alphaR, hs.alphaI, params.ssp->AttenUnit);
            hs.cS = crci(params, hsx.zTemp, hs.betaR, hs.betaI, params.ssp->AttenUnit);
            // printf("%g %g %g %g %c%c %g\n", zTemp, RecycledHS.alphaR,
            // RecycledHS.alphaI, params.freqinfo->freq0,
            //     params.ssp->AttenUnit[0], params.ssp->AttenUnit[1], params.fT);
            // printf("cp computed to (%g,%g)\n", hs.cP.real(), hs.cP.imag());
        } else if(hs.bc == 'G') { // *** Grain size (formulas from UW-APL HF Handbook)
            real vr, alpha2_f;    // values related to grain size

            // These formulas are from the UW-APL Handbook
            // The code is taken from older Matlab and is unnecesarily verbose
            // vr   is the sound speed ratio
            // rho is the density ratio

            if(hsx.Mz >= FL(-1.0) && hsx.Mz < FL(1.0)) {
                vr     = FL(0.002709) * SQ(hsx.Mz) - FL(0.056452) * hsx.Mz + FL(1.2778);
                hs.rho = FL(0.007797) * SQ(hsx.Mz) - FL(0.17057) * hsx.Mz + FL(2.3139);
            } else if(hsx.Mz >= FL(1.0) && hsx.Mz < FL(5.3)) {
                vr = FL(-0.0014881) * CUBE(hsx.Mz) + FL(0.0213937) * SQ(hsx.Mz)
                    - FL(0.1382798) * hsx.Mz + FL(1.3425);
                hs.rho = FL(-0.0165406) * CUBE(hsx.Mz) + FL(0.2290201) * SQ(hsx.Mz)
                    - FL(1.1069031) * hsx.Mz + FL(3.0455);
            } else {
                vr     = FL(-0.0024324) * hsx.Mz + FL(1.0019);
                hs.rho = FL(-0.0012973) * hsx.Mz + FL(1.1565);
            }

            if(hsx.Mz >= FL(-1.0) && hsx.Mz < FL(0.0)) {
                alpha2_f = FL(0.4556);
            } else if(hsx.Mz >= FL(0.0) && hsx.Mz < FL(2.6)) {
                alpha2_f = FL(0.4556) + FL(0.0245) * hsx.Mz;
            } else if(hsx.Mz >= FL(2.6) && hsx.Mz < FL(4.5)) {
                alpha2_f = FL(0.1978) + FL(0.1245) * hsx.Mz;
            } else if(hsx.Mz >= FL(4.5) && hsx.Mz < FL(6.0)) {
                alpha2_f = FL(8.0399) - FL(2.5228) * hsx.Mz + FL(0.20098) * SQ(hsx.Mz);
            } else if(hsx.Mz >= FL(6.0) && hsx.Mz < FL(9.5)) {
                alpha2_f = FL(0.9431) - FL(0.2041) * hsx.Mz + FL(0.0117) * SQ(hsx.Mz);
            } else {
                alpha2_f = FL(0.0601);
            }

            // params.ssp->AttenUnit = 'L';  // loss parameter
            // !! following uses a reference sound speed of 1500 ???
            // !! should be sound speed in the water, just above the sediment
            // the term vr / 1000 converts vr to units of m per ms
            hs.alphaR = vr * FL(1500.0);
            hs.alphaI = alpha2_f * (vr / FL(1000.0)) * FL(1500.0) * STD::log(FL(10.0))
                / (FL(40.0) * REAL_PI); // loss parameter Sect. IV., Eq. (4) of handbook

            hs.cP = crci(params, hsx.zTemp, hs.alphaR, hs.alphaI, {'L', ' '});
            hs.cS = FL(0.0);
        }
    }

private:
    BdryPtSmall &GetBdry(const bhcParams<O3D> &params) const
    {
        if constexpr(ISTOP)
            return params.Bdry->Top;
        else
            return params.Bdry->Bot;
    }
};

template<bool O3D> using BoundaryCondTop = BoundaryCond<O3D, true>;
template<bool O3D> using BoundaryCondBot = BoundaryCond<O3D, false>;

}} // namespace bhc::module
