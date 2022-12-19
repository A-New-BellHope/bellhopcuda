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
#pragma once
#include "common.hpp"

namespace bhc {

template<bool O3D> HOST_DEVICE inline bool IsRayRun(const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'R';
}

template<bool O3D> HOST_DEVICE inline bool IsTLRun(const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'C' || r == 'S' || r == 'I';
}

template<bool O3D> HOST_DEVICE inline bool IsEigenraysRun(const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'E';
}

template<bool O3D> HOST_DEVICE inline bool IsArrivalsRun(const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'A' || r == 'a';
}

template<bool O3D> HOST_DEVICE inline bool IsCoherentRun(const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'C';
}

template<bool O3D> HOST_DEVICE inline bool IsSemiCoherentRun(
    const BeamStructure<O3D> *Beam)
{
    char r = Beam->RunType[0];
    return r == 'S';
}

// Beam->Type[0] is
//   'G', '^', or ' ' Geometric hat beams in Cartesian coordinates
//   'g' Geometric hat beams in ray-centered coordinates
//   'B' Geometric Gaussian beams in Cartesian coordinates
//   'b' Geometric Gaussian beams in ray-centered coordinates
//   'S' Simple Gaussian beams
//   'C' Cerveny Gaussian beams in Cartesian coordinates
//   'R' Cerveny Gaussian beams in Ray-centered coordinates

template<bool O3D> HOST_DEVICE inline bool IsCervenyInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == 'R' || t == 'C';
}

template<bool O3D> HOST_DEVICE inline bool IsGeometricInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == ' ' || t == '^' || t == 'G' || t == 'g' || t == 'B' || t == 'b';
}

template<bool O3D> HOST_DEVICE inline bool IsSGBInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == 'S';
}

template<bool O3D> HOST_DEVICE inline bool IsCartesianInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == 'C' || t == ' ' || t == '^' || t == 'G' || t == 'B';
}

template<bool O3D> HOST_DEVICE inline bool IsRayCenInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == 'R' || t == 'g' || t == 'b';
}

template<bool O3D> HOST_DEVICE inline bool IsHatGeomInfl(const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == ' ' || t == '^' || t == 'G' || t == 'g';
}

template<bool O3D> HOST_DEVICE inline bool IsGaussianGeomInfl(
    const BeamStructure<O3D> *Beam)
{
    char t = Beam->Type[0];
    return t == 'B' || t == 'b';
}

template<bool O3D, bool R3D> HOST_DEVICE inline const char *GetBeamTypeTag(
    const BeamStructure<O3D> *Beam)
{
    switch(Beam->Type[0]) {
    case 'C':
    case 'R': return R3D ? "Cerveny style beam" : "Paraxial beams";
    case '^':
    case ' ':
        if constexpr(!R3D)
            return "Warning, Beam->Type[0] = ^ or ' ' not properly handled in BELLHOP "
                   "(2D)";
        [[fallthrough]];
    case 'G':
        if constexpr(R3D) return "Geometric beam, hat-shaped, Cart. coord.";
        [[fallthrough]];
    case 'g':
        if constexpr(R3D) return "Geometric beam, hat-shaped, Ray coord.";
        return "Geometric hat beams";
    case 'B':
        return R3D ? "Geometric beam, Gaussian-shaped, Cart. coord."
                   : "Geometric Gaussian beams";
    case 'b':
        return R3D
            ? "Geometric beam, Gaussian-shaped, Ray coord."
            : "Geo Gaussian beams in ray-cent. coords. not implemented in BELLHOP (2D)";
    case 'S': return "Simple Gaussian beams";
    default: GlobalLog("Invalid Beam->Type[0] %c\n", Beam->Type[0]); bail();
    }
}

// Beam->Type[1] controls the setting of the beam width
//   'F' space Filling
//   'M' minimum width
//   'W' WKB beams
//   LP: 'C': "Cerveny style beam"

template<bool O3D> HOST_DEVICE inline bool IsBeamWidthSpaceFilling(
    const BeamStructure<O3D> *Beam)
{
    return Beam->Type[1] == 'F';
}

template<bool O3D> HOST_DEVICE inline bool IsBeamWidthMinimum(
    const BeamStructure<O3D> *Beam)
{
    return Beam->Type[1] == 'M';
}

template<bool O3D> HOST_DEVICE inline bool IsBeamWidthWKB(const BeamStructure<O3D> *Beam)
{
    return Beam->Type[1] == 'W';
}

template<bool O3D> HOST_DEVICE inline bool IsBeamWidthCerveny(
    const BeamStructure<O3D> *Beam)
{
    return Beam->Type[1] == 'C';
}

template<bool O3D> HOST_DEVICE inline const char *GetBeamWidthTag(
    const BeamStructure<O3D> *Beam)
{
    switch(Beam->Type[1]) {
    case 'F': return "Space filling beams";
    case 'M': return "Minimum width beams";
    case 'W': return "WKB beams";
    case 'C': return "Cerveny style beam";
    default: GlobalLog("Invalid Beam->Type[1] %c\n", Beam->Type[1]); bail();
    }
}

// Beam->Type[2] controls curvature changes on boundary reflections
//   'D' Double
//   'S' Single
//   'Z' Zero
// Beam->Type[3] selects whether beam shifts are implemented on boundary reflection
//   'S' yes
//   'N' no

/**
 * LP: Defaults to 'R' (point source) for any other option in ReadRunType.
 */
template<bool O3D> HOST_DEVICE inline bool IsLineSource(const BeamStructure<O3D> *Beam)
{
    return Beam->RunType[3] == 'X';
}

/**
 * LP: Defaults to 'R' (rectilinear) for any other option in ReadRunType.
 */
template<bool O3D> HOST_DEVICE inline bool IsIrregularGrid(const BeamStructure<O3D> *Beam)
{
    return Beam->RunType[4] == 'I';
}

// ssp->Type:
//   'N': 1D   N2-linear profile option
//   'C': 1D   C-linear profile option
//   'S': 1D   Cubic spline profile option
//   'P': 1D   monotone PCHIP ACS profile option
//   'Q': 2D   Bilinear quadrilateral interpolation of SSP data in 2D
//   'H': 3D   Trilinear hexahedral interpolation of SSP data in 3D
//   'A': AnyD Analytic profile option

HOST_DEVICE inline bool IsN2LinearSSP(const SSPStructure *ssp)
{
    return ssp->Type == 'N';
}

HOST_DEVICE inline bool IsCLinearSSP(const SSPStructure *ssp) { return ssp->Type == 'C'; }

HOST_DEVICE inline bool IsCCubicSSP(const SSPStructure *ssp) { return ssp->Type == 'S'; }

HOST_DEVICE inline bool IsCPCHIPSSP(const SSPStructure *ssp) { return ssp->Type == 'P'; }

HOST_DEVICE inline bool IsQuadSSP(const SSPStructure *ssp) { return ssp->Type == 'Q'; }

HOST_DEVICE inline bool IsHexahedralSSP(const SSPStructure *ssp)
{
    return ssp->Type == 'H';
}

HOST_DEVICE inline bool IsAnalyticSSP(const SSPStructure *ssp)
{
    return ssp->Type == 'A';
}

HOST_DEVICE inline bool Is1DSSP(const SSPStructure *ssp)
{
    char t = ssp->Type;
    return t == 'N' || t == 'C' || t == 'S' || t == 'P';
}

HOST_DEVICE inline bool Is2DSSP(const SSPStructure *ssp)
{
    char t = ssp->Type;
    return t == 'Q';
}

HOST_DEVICE inline bool Is3DSSP(const SSPStructure *ssp)
{
    char t = ssp->Type;
    return t == 'H';
}

HOST_DEVICE inline bool IsAnyDSSP(const SSPStructure *ssp)
{
    char t = ssp->Type;
    return t == 'A';
}

} // namespace bhc
