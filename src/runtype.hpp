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

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common.hpp!"
#endif

namespace bhc {

template<char RT> struct RunType {
    static constexpr bool IsRay() { return RT == 'R'; }
    static constexpr bool IsTL() { return RT == 'C' /*|| RT == 'S' || RT == 'I'*/; }
    static constexpr bool IsEigenrays() { return RT == 'E'; }
    static constexpr bool IsArrivals() { return RT == 'A' /*|| RT == 'a'*/; }
    /*
    static constexpr bool IsCoherent() { return RT == 'C'; }
    static constexpr bool IsSemiCoherent() { return RT == 'S'; }
    static constexpr bool IsIncoherent() { return RT == 'I'; }
    */

    static_assert(
        IsRay() || IsTL() || IsEigenrays() || IsArrivals(),
        "RunType templated with invalid character!");
};

template<char IT> struct InflType {
    // ' ' and '^' are equivalent to 'G', but are handled in the template
    // selection, not here.
    static constexpr bool IsCerveny() { return IT == 'R' || IT == 'C'; }
    static constexpr bool IsGeometric() { return IT == 'G' || IT == 'g'; }
    // return IsHatGeom() || IsGaussianGeom();
    static constexpr bool IsSGB() { return IT == 'S'; }
    static constexpr bool IsCartesian()
    {
        return IT == 'C' || IT == 'G' /*|| IT == 'B'*/;
    }
    static constexpr bool IsRayCen() { return IT == 'R' || IT == 'g' /*|| IT == 'b'*/; }
    /*
    static constexpr bool IsHatGeom() { return IT == 'G' || IT == 'g'; }
    static constexpr bool IsGaussianGeom() { return IT == 'B' || IT == 'b'; }
    */

    static_assert(
        IsCerveny() || IsGeometric() || IsSGB(),
        "InflType templated with invalid character!");
};

template<char ST> struct SSPType {
    static constexpr bool IsN2Linear() { return ST == 'N'; }
    static constexpr bool IsCLinear() { return ST == 'C'; }
    static constexpr bool IsCCubic() { return ST == 'S'; }
    static constexpr bool IsCPCHIP() { return ST == 'P'; }
    static constexpr bool IsQuad() { return ST == 'Q'; }
    static constexpr bool IsHexahedral() { return ST == 'H'; }
    static constexpr bool IsAnalytic() { return ST == 'A'; }
    static constexpr bool Is1D()
    {
        return IsN2Linear() || IsCLinear() || IsCCubic() || IsCPCHIP();
    }
    static constexpr bool Is2D() { return IsQuad(); }
    static constexpr bool Is3D() { return IsHexahedral(); }
    static constexpr bool IsAnyD() { return IsAnalytic(); }

    static_assert(
        Is1D() || Is2D() || Is3D() || IsAnyD(),
        "SSPType templated with invalid character!");
};

template<char RT, char IT, char ST> struct CfgSel {
    using run  = RunType<RT>;
    using infl = InflType<IT>;
    using ssp  = SSPType<ST>;
};

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

template<bool O3D> HOST_DEVICE inline const char *GetBeamTypeTag(
    const BeamStructure<O3D> *Beam)
{
    switch(Beam->Type[0]) {
    case 'C':
    case 'R':
        if constexpr(O3D)
            return "Cerveny style beam";
        else
            return "Paraxial beams";
    case '^':
    case ' ':
        if constexpr(!O3D)
            return "Warning, Beam->Type[0] = ^ or ' ' not properly handled in BELLHOP "
                   "(2D)";
        [[fallthrough]];
    case 'G':
        if constexpr(O3D) return "Geometric beam, hat-shaped, Cart. coord.";
        [[fallthrough]];
    case 'g':
        if constexpr(O3D)
            return "Geometric beam, hat-shaped, Ray coord.";
        else
            return "Geometric hat beams";
    case 'B':
        if constexpr(O3D)
            return "Geometric beam, Gaussian-shaped, Cart. coord.";
        else
            return "Geometric Gaussian beams";
    case 'b':
        if constexpr(O3D)
            return "Geometric beam, Gaussian-shaped, Ray coord.";
        else
            return "Geo Gaussian beams in ray-cent. coords. not implemented in BELLHOP "
                   "(2D)";
    case 'S': return "Simple Gaussian beams";
    default: return "Invalid Beam->Type[0]";
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
    default: return "Invalid Beam->Type[1]";
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

} // namespace bhc
