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
    DefaultRcvrRanges<O3D, R3D>(params);
    DefaultRcvrBearings<O3D, R3D>(params);
    DefaultfreqVec<O3D, R3D>(params);
    DefaultRunType<O3D, R3D>(params);
    DefaultRayAnglesElevation<O3D, R3D>(params);
    DefaultRayAnglesBearing<O3D, R3D>(params);
    DefaultBeamInfo<O3D, R3D>(params);
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
    ReadRcvrRanges<O3D, R3D>(params, ENVFile);
    ReadRcvrBearings<O3D, R3D>(params, ENVFile);
    ReadfreqVec<O3D, R3D>(params, ENVFile);
    ReadRunType<O3D, R3D>(params, ENVFile);
    ReadRayAnglesElevation<O3D, R3D>(params, ENVFile);
    ReadRayAnglesBearing<O3D, R3D>(params, ENVFile);
    ReadBeamInfo<O3D, R3D>(params, ENVFile);

    ReadBoundary<O3D, R3D>(
        params, params.Bdry->Top.hs.Opt[4], params.Bdry->Top.hs.Depth,
        &params.bdinfo->top, true); // AlTImetry
    ReadBoundary<O3D, R3D>(
        params, params.Bdry->Bot.hs.Opt[1], params.Bdry->Bot.hs.Depth,
        &params.bdinfo->bot, false);             // BaThYmetry
    ReadReflectionCoefficient<O3D, R3D>(params); // (top and bottom)
    ReadPat<O3D, R3D>(params);                   // Source Beam Pattern

    GetInternal(params)->PRTFile << "\n";
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
