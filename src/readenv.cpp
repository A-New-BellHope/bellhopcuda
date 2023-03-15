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
#include "readenv.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

namespace bhc {

template<bool O3D, bool R3D> inline void PreprocessTopOpt(bhcParams<O3D, R3D> &params)
{
    params.ssp->Type         = params.Bdry->Top.hs.Opt[0];
    params.Bdry->Top.hs.bc   = params.Bdry->Top.hs.Opt[1];
    params.ssp->AttenUnit[0] = params.Bdry->Top.hs.Opt[2];
    params.ssp->AttenUnit[1] = params.Bdry->Top.hs.Opt[3];
}

template<bool O3D, bool R3D> inline void DefaultTopOpt(bhcParams<O3D, R3D> &params)
{
    // SSP (clinear), top bc (vacuum), atten units (dB/wavelength),
    // add vol atten (none), altimetry (none), dev mode (off)
    memcpy(params.Bdry->Top.hs.Opt, "CVW - ", 6);
    PreprocessTopOpt<O3D, R3D>(params);
}

/**
 * LP: Read top halfspace options; 4 out of the 6 entries are general program
 * options.
 */
template<bool O3D, bool R3D> inline void ReadTopOpt(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    LIST(ENVFile);
    ENVFile.Read(params.Bdry->Top.hs.Opt, 6); // LP: LDIFile fills rest with ' '
    PRTFile << "\n";

    PreprocessTopOpt<O3D, R3D>(params);

    // SSP approximation options

    switch(params.ssp->Type) {
    case 'N': PRTFile << "    N2-linear approximation to SSP\n"; break;
    case 'C': PRTFile << "    C-linear approximation to SSP\n"; break;
    case 'P': PRTFile << "    PCHIP approximation to SSP\n"; break;
    case 'S': PRTFile << "    Spline approximation to SSP\n"; break;
    case 'Q': {
        PRTFile << "    Quad approximation to SSP\n";
        // LP: This just checks for existence, moved actual open for reading
        // to InitQuad.
        std::ifstream SSPFile;
        SSPFile.open(GetInternal(params)->FileRoot + ".ssp");
        if(!SSPFile.good()) {
            PRTFile << "SSPFile = " << GetInternal(params)->FileRoot << ".ssp\n";
            EXTERR(BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file");
        }
    } break;
    case 'H': {
        // LP: This just checks for existence, moved actual open for reading
        // to InitHexahedral.
        PRTFile << "    Hexahedral approximation to SSP\n";
        std::ifstream SSPFile;
        SSPFile.open(GetInternal(params)->FileRoot + ".ssp");
        if(!SSPFile.good()) {
            PRTFile << "SSPFile = " << GetInternal(params)->FileRoot << ".ssp\n";
            EXTERR(BHC_PROGRAMNAME " - ReadEnvironment: Unable to open the SSP file");
        }
    } break;
    case 'A': PRTFile << "    Analytic SSP option\n"; break;
    default: EXTERR("ReadEnvironment: Unknown option for SSP approximation");
    }

    // Attenuation options

    switch(params.ssp->AttenUnit[0]) {
    case 'N': PRTFile << "    Attenuation units: nepers/m\n"; break;
    case 'F': PRTFile << "    Attenuation units: dB/mkHz\n"; break;
    case 'M': PRTFile << "    Attenuation units: dB/m\n"; break;
    case 'W': PRTFile << "    Attenuation units: dB/wavelength\n"; break;
    case 'Q': PRTFile << "    Attenuation units: Q\n"; break;
    case 'L': PRTFile << "    Attenuation units: Loss parameter\n"; break;
    default: EXTERR("ReadEnvironment: Unknown attenuation units");
    }

    // optional addition of volume attenuation using standard formulas

    AttenInfo *atten = params.atten;
    switch(params.ssp->AttenUnit[1]) {
    case 'T': PRTFile << "    THORP volume attenuation added\n"; break;
    case 'F':
        PRTFile << "    Francois-Garrison volume attenuation added\n";
        LIST(ENVFile);
        ENVFile.Read(atten->t);
        ENVFile.Read(atten->Salinity);
        ENVFile.Read(atten->pH);
        ENVFile.Read(atten->z_bar);
        PRTFile << std::setprecision(4);
        PRTFile << " T = " << std::setw(11) << atten->t
                << "degrees   S = " << std::setw(11) << atten->Salinity
                << " psu   pH = " << std::setw(11) << atten->pH
                << " z_bar = " << std::setw(11) << " m\n";
        break;
    case 'B':
        PRTFile << "    Biological attenuation\n";
        LIST(ENVFile);
        ENVFile.Read(atten->NBioLayers);
        PRTFile << "      Number of Bio Layers = " << atten->NBioLayers << "\n";
        for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio) {
            LIST(ENVFile);
            ENVFile.Read(atten->bio[iBio].z1);
            ENVFile.Read(atten->bio[iBio].z2);
            ENVFile.Read(atten->bio[iBio].f0);
            ENVFile.Read(atten->bio[iBio].q);
            ENVFile.Read(atten->bio[iBio].a0);
            PRTFile << "      Top    of layer = " << atten->bio[iBio].z1 << " m\n";
            PRTFile << "      Bottom of layer = " << atten->bio[iBio].z2 << " m\n";
            PRTFile << "      Resonance frequency = " << atten->bio[iBio].f0 << " Hz\n";
            PRTFile << "      Q = " << atten->bio[iBio].q << "\n";
            PRTFile << "      a0 = " << atten->bio[iBio].a0 << "\n";
        }
    case ' ': break;
    default: EXTERR("ReadEnvironment: Unknown top option letter in fourth position");
    }

    switch(params.Bdry->Top.hs.Opt[4]) {
    case '~':
    case '*': PRTFile << "    Altimetry file selected\n"; break;
    case '-':
    case '_':
    case ' ': break;
    default:
        PRTFile << "ReadEnvironment: Unknown top option letter in fifth position\n";
        std::abort();
    }

    switch(params.Bdry->Top.hs.Opt[5]) {
    case 'I': PRTFile << "    Development options enabled\n"; break;
    case ' ': break;
    default:
        PRTFile << "ReadEnvironment: Unknown top option letter in sixth position\n";
        std::abort();
    }
}

template<bool O3D, bool R3D> inline void PreprocessTopOpt(bhcParams<O3D, R3D> &params)
{
    params.Bdry->Bot.hs.bc = params.Bdry->Bot.hs.Opt[0];
}

template<bool O3D, bool R3D> inline void DefaultBotOpt(bhcParams<O3D, R3D> &params)
{
    memcpy(params.Bdry->Bot.hs.Opt, "R-    ", 6);
    PreprocessBotOpt<O3D, R3D>(params);
}

template<bool O3D, bool R3D> inline void ReadBotOpt(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    real Sigma;
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    LIST(ENVFile);
    ENVFile.Read(params.Bdry->Bot.hs.Opt, 6); // LP: LDIFile fills rest with ' '
    ENVFile.Read(Sigma);
    PRTFile << "\n RMS roughness = " << std::setw(10) << std::setprecision(3) << Sigma
            << "\n";

    switch(params.Bdry->Bot.hs.Opt[1]) {
    case '~':
    case '*': PRTFile << "    Bathymetry file selected\n"; break;
    case '-':
    case '_':
    case ' ': break;
    default:
        EXTERR(
            "Unknown bottom option letter in second position: Bdry->Bot.hs.Opt[1] == "
            "'%c'",
            params.Bdry->Bot.hs.Opt[1]);
    }

    PreprocessBotOpt<O3D, R3D>(params);
}

template<bool O3D, bool R3D> inline void DefaultRunType(bhcParams<O3D, R3D> &params)
{
    // RunType, infl/beam type, ignored, point source, rectilinear grid, dim, ignored
    memcpy(params.Beam->RunType, R3D ? "CG RR3 " : "CG RR2 ", 7);
}

/**
 * Read the RunType variable and echo with explanatory information to the print file
 */
template<bool O3D, bool R3D> inline void ReadRunType(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    LIST(ENVFile);
    ENVFile.Read(params.Beam->RunType, 7);
    PRTFile << "\n";

    switch(params.Beam->RunType[0]) {
    case 'R': PRTFile << "Ray trace run\n"; break;
    case 'E': PRTFile << "Eigenray trace run\n"; break;
    case 'I': PRTFile << "Incoherent TL calculation\n"; break;
    case 'S': PRTFile << "Semi-coherent TL calculation\n"; break;
    case 'C': PRTFile << "Coherent TL calculation\n"; break;
    case 'A': PRTFile << "Arrivals calculation, ASCII  file output\n"; break;
    case 'a': PRTFile << "Arrivals calculation, binary file output\n"; break;
    default: EXTERR("ReadEnvironment: Unknown RunType selected");
    }

    switch(params.Beam->RunType[1]) {
    case 'C': PRTFile << "Cartesian beams\n"; break;
    case 'R': PRTFile << "Ray centered beams\n"; break;
    case 'S': PRTFile << "Simple gaussian beams\n"; break;
    case 'b': PRTFile << "Geometric gaussian beams in ray-centered coordinates\n"; break;
    case 'B': PRTFile << "Geometric gaussian beams in Cartesian coordinates\n"; break;
    case 'g': PRTFile << "Geometric hat beams in ray-centered coordinates\n"; break;
    default:
        params.Beam->RunType[1] = 'G';
        PRTFile << "Geometric hat beams in Cartesian coordinates\n";
    }

    switch(params.Beam->RunType[3]) {
    case 'X': PRTFile << "Line source (Cartesian coordinates)\n"; break;
    default:
        params.Beam->RunType[3] = 'R';
        PRTFile << "Point source (cylindrical coordinates)\n";
    }

    switch(params.Beam->RunType[4]) {
    case 'I':
        PRTFile << "Irregular grid: Receivers at Rr[:] x Rz[:]\n";
        if(params.Pos->NRz != params.Pos->NRr)
            EXTWARN("ReadEnvironment: Irregular grid option selected with NRz not "
                    "equal to Nr");
        // memcpy(PlotType, "irregular ", 10);
        break;
    default:
        PRTFile << "Rectilinear receiver grid: Receivers at Rr[:] x Rz[:]\n";
        params.Beam->RunType[4] = 'R';
        // memcpy(PlotType, "rectilin  ", 10);
    }

    bool defaulted2 = false;
    if(params.Beam->RunType[5] != '2' && params.Beam->RunType[5] != '3') {
        if constexpr(R3D) {
            EXTERR("Environment file does not specify dimensionality, defaults to 2 "
                   "(2D or Nx2D), but you are running " BHC_PROGRAMNAME " in 3D mode");
        } else if constexpr(O3D) {
            EXTWARN("Environment file does not specify dimensionality, defaults to 2 "
                    "(2D or Nx2D), assuming this is Nx2D because that is what you're "
                    "running");
        }
        params.Beam->RunType[5] = '2';
        defaulted2              = true;
    }
    if(params.Beam->RunType[5] == '2') {
        if(!defaulted2) {
            PRTFile << "N x 2D calculation (neglects horizontal refraction)\n";
        }
        if constexpr(R3D) {
            EXTERR("This is a 2D or Nx2D environment file, but you are "
                   "running " BHC_PROGRAMNAME " in 3D mode");
        }
    } else {
        PRTFile << "3D calculation\n";
        if constexpr(!R3D) {
            EXTERR("This is a 3D environment file, but you are running " BHC_PROGRAMNAME
                   " in 2D or Nx2D mode");
        }
    }
}

template<bool O3D, bool R3D> inline void DefaultBoundaryCond(
    bhcParams<O3D, R3D> &params, HSInfo &hs)
{
    hs.cP = hs.cS = hs.rho = FL(0.0);
}

/**
 * LP: Formerly TopBot
 * Handles top and bottom boundary conditions
 * params.freqinfo->freq0: center / nominal frequency (wideband not supported)
 */
template<bool O3D, bool R3D> inline void ReadBoundaryCond(
    bhcParams<O3D, R3D> &params, HSInfo &hs, LDIFile &ENVFile, HSInfo &RecycledHS,
    bool isTop)
{
    real Mz, vr, alpha2_f; // values related to grain size
    real zTemp;

    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    // LP: Moved from ReadEnvironment.
    if(isTop && hs.bc == 'A') {
        PRTFile << "      z         alphaR      betaR     rho        alphaI     betaI\n";
        PRTFile << "     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)\n";
    }

    // Echo to PRTFile user's choice of boundary condition

    switch(hs.bc) {
    case 'V': PRTFile << "    VACUUM\n"; break;
    case 'R': PRTFile << "    Perfectly RIGID\n"; break;
    case 'A': PRTFile << "    ACOUSTO-ELASTIC half-space\n"; break;
    case 'G': PRTFile << "    Grain size to define half-space\n"; break;
    case 'F': PRTFile << "    FILE used for reflection loss\n"; break;
    case 'W': PRTFile << "    Writing an IFL file\n"; break;
    case 'P': PRTFile << "    reading PRECALCULATED IFL\n"; break;
    default: EXTERR("TopBot: Unknown boundary condition type");
    }

    // ****** Read in BC parameters depending on particular choice ******

    hs.cP = hs.cS = hs.rho = FL(0.0);

    if(hs.bc == 'A') { // *** Half-space properties ***
        zTemp = FL(0.0);
        LIST(ENVFile);
        ENVFile.Read(zTemp);
        ENVFile.Read(RecycledHS.alphaR);
        ENVFile.Read(RecycledHS.betaR);
        ENVFile.Read(RecycledHS.rho);
        ENVFile.Read(RecycledHS.alphaI);
        ENVFile.Read(RecycledHS.betaI);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " " << std::setw(10)
                << RecycledHS.alphaR << " " << std::setw(10) << RecycledHS.betaR << " "
                << std::setw(6) << RecycledHS.rho << " " << std::setprecision(4)
                << std::setw(10) << RecycledHS.alphaI << " " << std::setw(10)
                << RecycledHS.betaI << "\n";
        // dummy parameters for a layer with a general power law for attenuation
        // these are not in play because the AttenUnit for this is not allowed yet
        // freq0         = freq;
        params.fT = FL(1000.0);

        hs.cP = crci(
            params, zTemp, RecycledHS.alphaR, RecycledHS.alphaI, params.ssp->AttenUnit);
        hs.cS = crci(
            params, zTemp, RecycledHS.betaR, RecycledHS.betaI, params.ssp->AttenUnit);
        // printf("%g %g %g %g %c%c %g\n", zTemp, RecycledHS.alphaR,
        // RecycledHS.alphaI, params.freqinfo->freq0,
        //     params.ssp->AttenUnit[0], params.ssp->AttenUnit[1], params.fT);
        // printf("cp computed to (%g,%g)\n", hs.cP.real(), hs.cP.imag());

        hs.rho = RecycledHS.rho;
    } else if(hs.bc == 'G') { // *** Grain size (formulas from UW-APL HF Handbook)

        // These formulas are from the UW-APL Handbook
        // The code is taken from older Matlab and is unnecesarily verbose
        // vr   is the sound speed ratio
        // rho is the density ratio
        LIST(ENVFile);
        ENVFile.Read(zTemp);
        ENVFile.Read(Mz);
        PRTFile << std::setprecision(2) << std::setw(10) << zTemp << " " << std::setw(10)
                << Mz << "\n";

        if(Mz >= FL(-1.0) && Mz < FL(1.0)) {
            vr             = FL(0.002709) * SQ(Mz) - FL(0.056452) * Mz + FL(1.2778);
            RecycledHS.rho = FL(0.007797) * SQ(Mz) - FL(0.17057) * Mz + FL(2.3139);
        } else if(Mz >= FL(1.0) && Mz < FL(5.3)) {
            vr = FL(-0.0014881) * CUBE(Mz) + FL(0.0213937) * SQ(Mz) - FL(0.1382798) * Mz
                + FL(1.3425);
            RecycledHS.rho = FL(-0.0165406) * CUBE(Mz) + FL(0.2290201) * SQ(Mz)
                - FL(1.1069031) * Mz + FL(3.0455);
        } else {
            vr             = FL(-0.0024324) * Mz + FL(1.0019);
            RecycledHS.rho = FL(-0.0012973) * Mz + FL(1.1565);
        }

        if(Mz >= FL(-1.0) && Mz < FL(0.0)) {
            alpha2_f = FL(0.4556);
        } else if(Mz >= FL(0.0) && Mz < FL(2.6)) {
            alpha2_f = FL(0.4556) + FL(0.0245) * Mz;
        } else if(Mz >= FL(2.6) && Mz < FL(4.5)) {
            alpha2_f = FL(0.1978) + FL(0.1245) * Mz;
        } else if(Mz >= FL(4.5) && Mz < FL(6.0)) {
            alpha2_f = FL(8.0399) - FL(2.5228) * Mz + FL(0.20098) * SQ(Mz);
        } else if(Mz >= FL(6.0) && Mz < FL(9.5)) {
            alpha2_f = FL(0.9431) - FL(0.2041) * Mz + FL(0.0117) * SQ(Mz);
        } else {
            alpha2_f = FL(0.0601);
        }

        // params.ssp->AttenUnit = 'L';  // loss parameter
        // !! following uses a reference sound speed of 1500 ???
        // !! should be sound speed in the water, just above the sediment
        // the term vr / 1000 converts vr to units of m per ms
        RecycledHS.alphaR = vr * FL(1500.0);
        RecycledHS.alphaI = alpha2_f * (vr / FL(1000.0)) * FL(1500.0) * STD::log(FL(10.0))
            / (FL(40.0) * REAL_PI); // loss parameter Sect. IV., Eq. (4) of handbook

        hs.cP  = crci(params, zTemp, RecycledHS.alphaR, RecycledHS.alphaI, {'L', ' '});
        hs.cS  = FL(0.0);
        hs.rho = RecycledHS.rho;
    }
}

template<bool O3D, bool R3D> inline void SetTitle(
    bhcParams<O3D, R3D> &params, const std::string &TempTitle)
{
    size_t l = bhc::min(sizeof(params.Title) - 1, TempTitle.size());
    memcpy(params.Title, TempTitle.c_str(), l);
    params.Title[l] = 0;
}

template<bool O3D, bool R3D> inline void ReadTitle(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    // Prepend model name to title
    std::string TempTitle;
    LIST(ENVFile);
    ENVFile.Read(TempTitle);
    TempTitle = BHC_PROGRAMNAME "- " + TempTitle;
    PRTFile << TempTitle << "\n";
    SetTitle<O3D, R3D>(params, TempTitle);
}

template<bool O3D, bool R3D> inline void DefaultTitle(bhcParams<O3D, R3D> &params)
{
    SetTitle<O3D, R3D>(params, BHC_PROGRAMNAME "- no env file");
}

template<bool O3D, bool R3D> inline void ReadFreq0(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    LIST(ENVFile);
    ENVFile.Read(params.freqinfo->freq0);
    PRTFile << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    PRTFile << " frequency = " << std::setw(11) << params.freqinfo->freq0 << " Hz\n";
}

template<bool O3D, bool R3D> inline void DefaultFreq0(bhcParams<O3D, R3D> &params)
{
    params.freqinfo->freq0 = RL(50.0);
}

template<bool O3D, bool R3D> inline void ReadNMedia(
    bhcParams<O3D, R3D> &params, LDIFile &ENVFile)
{
    int32_t NMedia;
    LIST(ENVFile);
    ENVFile.Read(NMedia);
    PRTFile << "Dummy parameter NMedia = " << NMedia << "\n";
    if(NMedia != 1) {
        EXTWARN("ReadEnvironment: Only one medium or layer is allowed in BELLHOP; "
                "sediment layers must be handled using a reflection coefficient");
    }
}

template<bool O3D, bool R3D> inline void TopDepthFromSSP(bhcParams<O3D, R3D> &params)
{
    // Depth of top boundary is taken from first SSP point
    // LP: Must be here, used later in env setup.
    Bdry->Top.hs.Depth = params.ssp->z[0];
    // [mbp:] bottom depth should perhaps be set the same way?
    // TODO some preprocessing to check that SSP bounds equal Top/Bot bounds
}

template<bool O3D, bool R3D> void DefaultEnvironment(bhcParams<O3D, R3D> &params)
{
    DefaultTitle<O3D, R3D>(params);
    DefaultFreq0<O3D, R3D>(params);
    DefaultTopOpt<O3D, R3D>(params);
    DefaultBoundaryCond<O3D, R3D>(params, params.Bdry->Top.hs);
    DefaultSSP<O3D, R3D>(params);
    TopDepthFromSSP<O3D, R3D>(params);
    DefaultBotOpt<O3D, R3D>(params);
    DefaultBoundaryCond<O3D, R3D>(params, params.Bdry->Bot.hs);
    DefaultSxSy<O3D, R3D>(params);
    DefaultSzRz<O3D, R3D>(params);
    DefaultRcvrRanges(params, ENVFile);
    if constexpr(O3D) DefaultRcvrBearings(params, ENVFile);
}

template<bool O3D, bool R3D> void ReadEnvironment(bhcParams<O3D, R3D> &params)
{
    PrintFileEmu &PRTFile = GetInternal(params)->PRTFile;

    // Values only initialized once--reused from top to ssp, and ssp to bot
    HSInfo RecycledHS;
    RecycledHS.alphaR = FL(1500.0);
    RecycledHS.betaR  = FL(0.0);
    RecycledHS.alphaI = FL(0.0);
    RecycledHS.betaI  = FL(0.0);
    RecycledHS.rho    = FL(1.0);

    PRTFile << BHC_PROGRAMNAME << (R3D ? "3D" : O3D ? "Nx2D" : "") << "\n\n";

    // Open the environmental file
    LDIFile ENVFile(GetInternal(params), GetInternal(params)->FileRoot + ".env");
    if(!ENVFile.Good()) {
        PRTFile << "ENVFile = " << GetInternal(params)->FileRoot << ".env\n";
        EXTERR(BHC_PROGRAMNAME
               " - ReadEnvironment: Unable to open the environmental file");
    }

    ReadTitle<O3D, R3D>(params, ENVFile);
    ReadFreq0<O3D, R3D>(params, ENVFile);
    ReadNMedia<O3D, R3D>(params, ENVFile);
    ReadTopOpt<O3D, R3D>(params, ENVFile);
    ReadBoundaryCond<O3D, R3D>(params, params.Bdry->Top.hs, ENVFile, RecycledHS, true);
    ReadSSP<O3D, R3D>(params, ENVFile, RecycledHS);
    TopDepthFromSSP<O3D, R3D>(params);
    ReadBotOpt<O3D, R3D>(params, ENVFile);
    ReadBoundaryCond<O3D, R3D>(params, params.Bdry->Bot.hs, ENVFile, RecycledHS, false);
    ReadSxSy<O3D, R3D>(params, ENVFile);
    ReadSzRz<O3D, R3D>(params, ENVFile);
    ReadRcvrRanges(params, ENVFile);
    if constexpr(O3D) ReadRcvrBearings(params, ENVFile);

    ReadfreqVec(params, ENVFile);
    ReadRunType<O3D, R3D>(params, ENVFile);

    real Depth = params.Bdry->Bot.hs.Depth - params.Bdry->Top.hs.Depth; // water depth
    ReadRayAngles<O3D, R3D, false>(params, Depth, ENVFile, params.Angles->alpha);
    if constexpr(O3D) {
        ReadRayAngles<O3D, R3D, true>(params, Depth, ENVFile, params.Angles->beta);
    }

    PRTFile << "\n_______________________________________________________________________"
               "___\n\n";

    // LP: Moved to separate function for clarity and modularity.
    ReadBeamInfo<O3D, R3D>(params, ENVFile);
}

#if BHC_ENABLE_2D
template void ReadEnvironment<false, false>(
    bhcParams<false, false> &params, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_NX2D
template void ReadEnvironment<true, false>(
    bhcParams<true, false> &params, HSInfo &RecycledHS);
#endif
#if BHC_ENABLE_3D
template void ReadEnvironment<true, true>(
    bhcParams<true, true> &params, HSInfo &RecycledHS);
#endif

} // namespace bhc
