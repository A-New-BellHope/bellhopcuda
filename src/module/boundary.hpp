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
#include "../boundary.hpp"

namespace bhc { namespace module {

/**
 * Templated to become Altimetry or Bathymetry
 */
template<bool O3D, bool ISTOP> class Boundary : public ParamsModule<O3D> {
public:
    Boundary() {}
    virtual ~Boundary() {}

    virtual void Init(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        bdinfotb->bd                  = nullptr;
    }

    virtual void SetupPre(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        if constexpr(O3D) {
            bdinfotb->NPts.x = 2;
            bdinfotb->NPts.y = 2;
        } else {
            bdinfotb->NPts = 2;
        }
        memcpy(bdinfotb->type, "LS", 2);
        bdinfotb->rangeInKm = true;
        bdinfotb->dirty     = true;
    }

    virtual void Default(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);

        if constexpr(O3D) {
            bdinfotb->type[0] = 'R';
            bdinfotb->NPts    = int2(2, 2);
            trackallocate(params, s_altimetrybathymetry, bdinfotb->bd, 2 * 2);

            bdinfotb->bd[0].x = vec3(-BDRYBIG, -BDRYBIG, BdryDepth(params));
            bdinfotb->bd[1].x = vec3(-BDRYBIG, BDRYBIG, BdryDepth(params));
            bdinfotb->bd[2].x = vec3(BDRYBIG, -BDRYBIG, BdryDepth(params));
            bdinfotb->bd[3].x = vec3(BDRYBIG, BDRYBIG, BdryDepth(params));

            for(int32_t i = 0; i < 4; ++i) {
                bdinfotb->bd[i].t  = vec3(FL(1.0), FL(0.0), FL(0.0));
                bdinfotb->bd[i].n1 = vec3(FL(0.0), FL(0.0), NegTop);
                bdinfotb->bd[i].n2 = vec3(FL(0.0), FL(0.0), NegTop);
            }
            bdinfotb->dirty = false; // LP: No ComputeBdryTangentNormal cause done
                                     // manually here
        } else {
            trackallocate(params, s_altimetrybathymetry, bdinfotb->bd, 2);
            bdinfotb->bd[0].x = vec2(-BDRYBIG, BdryDepth(params));
            bdinfotb->bd[1].x = vec2(BDRYBIG, BdryDepth(params));
        }
    }

    virtual void Read(bhcParams<O3D> &params, LDIFile &, HSInfo &) const override
    {
        if(!IsFile(params)) {
            Default(params);
            return;
        }

        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        PrintFileEmu &PRTFile         = GetInternal(params)->PRTFile;

        LDIFile
            BDRYFile(GetInternal(params), GetInternal(params)->FileRoot + "." + s_atibty);
        if(!BDRYFile.Good()) {
            PRTFile << s_ATIBTY << "File = " << GetInternal(params)->FileRoot << "."
                    << s_atibty << "\n";
            EXTERR("Read%s: Unable to open %s file", s_ATIBTY, s_altimetrybathymetry);
        }

        LIST(BDRYFile);
        BDRYFile.Read(bdinfotb->type, O3D ? 1 : 2);
        SetupPost(params);

        if constexpr(O3D) {
            real *Globalx = nullptr, *Globaly = nullptr;

            // x values
            LIST(BDRYFile);
            BDRYFile.Read(bdinfotb->NPts.x);
            trackallocate(
                params, "temp arrays for altimetry/bathymetry", Globalx,
                bhc::max(bdinfotb->NPts.x, 3));
            Globalx[2] = FL(-999.9);
            LIST(BDRYFile);
            BDRYFile.Read(Globalx, bdinfotb->NPts.x);
            SubTab(Globalx, bdinfotb->NPts.x);

            // y values
            LIST(BDRYFile);
            BDRYFile.Read(bdinfotb->NPts.y);
            trackallocate(
                params, "temp arrays for altimetry/bathymetry", Globaly,
                bhc::max(bdinfotb->NPts.y, 3));
            Globaly[2] = FL(-999.9);
            LIST(BDRYFile);
            BDRYFile.Read(Globaly, bdinfotb->NPts.y);
            SubTab(Globaly, bdinfotb->NPts.y);

            // z values
            trackallocate(
                params, s_altimetrybathymetry, bdinfotb->bd,
                bdinfotb->NPts.x * bdinfotb->NPts.y);

            for(int32_t iy = 0; iy < bdinfotb->NPts.y; ++iy) {
                LIST(BDRYFile); // read a row of depths
                for(int32_t ix = 0; ix < bdinfotb->NPts.x; ++ix) {
                    vec3 &x = bdinfotb->bd[ix * bdinfotb->NPts.y + iy].x;
                    BDRYFile.Read(x.z);
                    x.x = Globalx[ix];
                    x.y = Globaly[iy];
                }
            }

            trackdeallocate(params, Globalx);
            trackdeallocate(params, Globaly);
        } else {
            LIST(BDRYFile);
            BDRYFile.Read(bdinfotb->NPts);
            bdinfotb->NPts += 2; // we'll be extending the s_altimetrybathymetry to
                                 // infinity to the left and right

            trackallocate(params, s_altimetrybathymetry, bdinfotb->bd, bdinfotb->NPts);

            for(int32_t ii = 1; ii < bdinfotb->NPts - 1; ++ii) {
                LIST(BDRYFile);
                BDRYFile.Read(bdinfotb->bd[ii].x);
                if(bdinfotb->type[1] == 'L') {
                    BDRYFile.Read(bdinfotb->bd[ii].hs.alphaR);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.betaR);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.rho);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.alphaI);
                    BDRYFile.Read(bdinfotb->bd[ii].hs.betaI);
                }
            }

            // extend the bathymetry to +/- infinity in a piecewise constant fashion
            bdinfotb->bd[0].x[0]                  = -BDRYBIG * RL(0.001); // in km
            bdinfotb->bd[0].x[1]                  = bdinfotb->bd[1].x[1];
            bdinfotb->bd[0].hs                    = bdinfotb->bd[1].hs;
            bdinfotb->bd[bdinfotb->NPts - 1].x[0] = BDRYBIG * RL(0.001); // in km
            bdinfotb->bd[bdinfotb->NPts - 1].x[1] = bdinfotb->bd[bdinfotb->NPts - 2].x[1];
            bdinfotb->bd[bdinfotb->NPts - 1].hs   = bdinfotb->bd[bdinfotb->NPts - 2].hs;
        }
    }

    virtual void Write(bhcParams<O3D> &params, LDOFile &) const
    {
        if(!IsFile(params)) return;

        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        PrintFileEmu &PRTFile         = GetInternal(params)->PRTFile;

        LDOFile BDRYFile;
        int32_t NPts;
        if constexpr(O3D)
            NPts = bdinfotb->NPts.x * bdinfotb->NPts.y;
        else
            NPts = bdinfotb->NPts;
        BDRYFile.setStyle(
            NPts <= 10 ? LDOFile::Style::WRITTEN_BY_HAND : LDOFile::Style::MATLAB_OUTPUT);
        BDRYFile.open(GetInternal(params)->FileRoot + "." + s_atibty);
        if(!BDRYFile.good()) {
            PRTFile << s_ATIBTY << "File = " << GetInternal(params)->FileRoot << "."
                    << s_atibty << "\n";
            EXTERR("Read%s: Unable to open new %s file", s_ATIBTY, s_altimetrybathymetry);
        }

        BDRYFile << std::string(bdinfotb->type, 2);
        BDRYFile.write("! " BHC_PROGRAMNAME "- ");
        BDRYFile.write(s_altimetrybathymetry);
        BDRYFile.write(" file for ");
        BDRYFile.write(params.Title);
        BDRYFile << '\n';

        if constexpr(O3D) {
            // x values
            UnSubTab(
                BDRYFile, &bdinfotb->bd[0].x.x, bdinfotb->NPts.x, "NPts.x", nullptr,
                RL(0.001),
                bdinfotb->NPts.y * sizeof(bdinfotb->bd[0]) / sizeof(bdinfotb->bd[0].x.x));

            // y values
            UnSubTab(
                BDRYFile, &bdinfotb->bd[0].x.y, bdinfotb->NPts.y, "NPts.y", nullptr,
                RL(0.001), sizeof(bdinfotb->bd[0]) / sizeof(bdinfotb->bd[0].x.y));

            // z values
            for(int32_t iy = 0; iy < bdinfotb->NPts.y; ++iy) {
                for(int32_t ix = 0; ix < bdinfotb->NPts.x; ++ix) {
                    vec3 &x = bdinfotb->bd[ix * bdinfotb->NPts.y + iy].x;
                    BDRYFile << x.z;
                }
                BDRYFile << '\n';
            }
        } else {
            BDRYFile << bdinfotb->NPts - 2;
            BDRYFile.write("! NPts\n");

            for(int32_t ii = 1; ii < bdinfotb->NPts - 1; ++ii) {
                BDRYFile << bdinfotb->bd[ii].x.x * RL(0.001) << bdinfotb->bd[ii].x.y;
                if(bdinfotb->type[1] == 'L') {
                    BDRYFile << bdinfotb->bd[ii].hs.alphaR;
                    BDRYFile << bdinfotb->bd[ii].hs.betaR;
                    BDRYFile << bdinfotb->bd[ii].hs.rho;
                    BDRYFile << bdinfotb->bd[ii].hs.alphaI;
                    BDRYFile << bdinfotb->bd[ii].hs.betaI;
                }
                BDRYFile << '\n';
            }
        }
    }

    virtual void SetupPost(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        if constexpr(O3D) bdinfotb->type[1] = ' ';
    }

    void ExtSetup(bhcParams<O3D> &params, const IORI2<O3D> &NPts) const
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        GetModeFlag(params)           = '~';
        bdinfotb->dirty               = true;
        bdinfotb->NPts                = NPts;
        if constexpr(O3D) {
            trackallocate(
                params, s_altimetrybathymetry, bdinfotb->bd,
                bdinfotb->NPts.x * bdinfotb->NPts.y);
        } else {
            trackallocate(params, s_altimetrybathymetry, bdinfotb->bd, bdinfotb->NPts);
        }
    }

    virtual void Validate(bhcParams<O3D> &params) const override
    {
        PrintFileEmu &PRTFile         = GetInternal(params)->PRTFile;
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);

        switch(bdinfotb->type[0]) {
        case 'R':
            if constexpr(!O3D) {
                EXTERR("%sType R not supported for 2D runs\n", s_atibty);
            }
            break;
        case 'C': break;
        case 'L':
            if constexpr(O3D) {
                EXTERR("%sType L not supported for 3D runs\n", s_atibty);
            }
            break;
        default:
            EXTERR(
                "Read%s: Unknown option for selecting %s interpolation", s_ATIBTY,
                s_altimetrybathymetry);
        }

        if constexpr(O3D) {
            if(bdinfotb->type[1] != ' ') {
                EXTERR(
                    "Read%s: %s option (type[1]) must be ' ' in Nx2D/3D mode\n", s_ATIBTY,
                    s_altimetrybathymetry);
            }
            if(!monotonic(
                   &bdinfotb->bd[0].x.x, bdinfotb->NPts.x,
                   bdinfotb->NPts.y * BdryStride<O3D>, 0)) {
                EXTERR(
                    "BELLHOP:Read%s: %s x-coordinates are not monotonically increasing",
                    s_ATIBTY, s_AltimetryBathymetry);
            }
            if(!monotonic(&bdinfotb->bd[0].x.y, bdinfotb->NPts.y, BdryStride<O3D>, 0)) {
                EXTERR(
                    "BELLHOP:Read%s: %s y-coordinates are not monotonically increasing",
                    s_ATIBTY, s_AltimetryBathymetry);
            }

            bool warnedNaN = false;
            for(int32_t iy = 0; iy < bdinfotb->NPts.y; ++iy) {
                for(int32_t ix = 0; ix < bdinfotb->NPts.x; ++ix) {
                    const vec3 &x = bdinfotb->bd[ix * bdinfotb->NPts.y + iy].x;
                    if(!std::isfinite(x.z) && !warnedNaN) {
                        PRTFile << "Warning in " BHC_PROGRAMNAME "3D - Read" << s_ATIBTY
                                << "3D : The " << s_altimetrybathymetry
                                << " file contains a NaN\n";
                        warnedNaN = true;
                    }
                }
            }
        } else {
            if(bdinfotb->type[1] != 'S' && bdinfotb->type[1] != ' '
               && bdinfotb->type[1] != 'L') {
                EXTERR(
                    "Read%s: Unknown option for selecting %s option", s_ATIBTY,
                    s_altimetrybathymetry);
            }

            for(int32_t ii = 0; ii < bdinfotb->NPts; ++ii) {
                if(-NegTop * bdinfotb->bd[ii].x[1] < -NegTop * BdryDepth(params)) {
                    EXTERR(
                        "BELLHOP:Read%s: %s %s point in the sound speed profile",
                        s_ATIBTY, s_AltimetryBathymetry, s_risesdrops);
                }
            }

            if(!monotonic(&bdinfotb->bd[0].x.x, bdinfotb->NPts, BdryStride<O3D>, 0)) {
                EXTERR(
                    "BELLHOP:Read%s: %s ranges are not monotonically increasing\n",
                    s_ATIBTY, s_AltimetryBathymetry);
            }
        }
    }

    virtual void Echo(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        PrintFileEmu &PRTFile         = GetInternal(params)->PRTFile;

        Preprocess(params);

        if(IsFile(params)) {
            PRTFile
                << "_____________________________________________________________________"
                   "_____\n\n";
            PRTFile << "Using " << s_topbottom << "-" << s_altimetrybathymetry
                    << " file\n";
        }
        // LP: Normally the rest of this is not echoed if it's not a file, but
        // we might want to echo user-set values.

        switch(bdinfotb->type[0]) {
        case 'R': PRTFile << "Regular grid for a 3D run\n"; break;
        case 'C':
            if constexpr(O3D) {
                PRTFile << "Regular grid for a 3D run (curvilinear)\n";
            } else {
                PRTFile << "Curvilinear Interpolation\n";
            }
            break;
        case 'L': PRTFile << "Piecewise linear interpolation\n"; break;
        }

        if constexpr(O3D) {
            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in x "
                    << bdinfotb->NPts.x << "\n";
            EchoVector(
                &bdinfotb->bd[0].x.x, bdinfotb->NPts.x, PRTFile, Bdry_Number_to_Echo, "",
                RL(0.001), bdinfotb->NPts.y * BdryStride<O3D>, 0);

            PRTFile << "\nNumber of " << s_altimetrybathymetry << " points in y "
                    << bdinfotb->NPts.y << "\n";
            EchoVector(
                &bdinfotb->bd[0].x.y, bdinfotb->NPts.y, PRTFile, Bdry_Number_to_Echo, "",
                RL(0.001), BdryStride<O3D>, 0);

            PRTFile << "\n";
        } else {
            PRTFile << "Number of " << s_altimetrybathymetry
                    << " points = " << bdinfotb->NPts << "\n";

            // LP: BUG: Geoacoustics are supported for altimetry, but the
            // header for geoacoustics is only supported for bathymetry.
            bool shortFormat = true;
            if constexpr(!ISTOP) {
                shortFormat = bdinfotb->type[1] == 'S' || bdinfotb->type[1] == ' ';
            }
            if(shortFormat) {
                if constexpr(!ISTOP) {
                    PRTFile << "Short format (" << s_altimetrybathymetry << " only)\n";
                }
                PRTFile << "\n Range (km)  Depth (m)\n";
            } else if(bdinfotb->type[1] == 'L') {
                PRTFile << "Long format (" << s_altimetrybathymetry
                        << " and geoacoustics)\n";
                PRTFile << "Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  "
                           "alphaI     betaI\n";
            }

            for(int32_t ii = 1; ii < bdinfotb->NPts - 1; ++ii) {
                // LP: This condition was previously ii == bdinfotb->NPts - 1,
                // which will never be satisfied due to the loop bounds
                // echo some values
                if(ii < Bdry_Number_to_Echo || ii == bdinfotb->NPts - 2) {
                    vec2 x = bdinfotb->bd[ii].x;
                    x.x *= RL(0.001); // convert back to km
                    PRTFile << std::setprecision(3) << x;
                    if(bdinfotb->type[1] == 'L') {
                        PRTFile << " " << bdinfotb->bd[ii].hs.alphaR << " "
                                << bdinfotb->bd[ii].hs.betaR << " "
                                << bdinfotb->bd[ii].hs.rho << " "
                                << bdinfotb->bd[ii].hs.alphaI << " "
                                << bdinfotb->bd[ii].hs.betaI;
                    }
                    PRTFile << "\n";
                }
            }
        }
    }

    virtual void Preprocess(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);

        if(bdinfotb->rangeInKm) {
            bdinfotb->rangeInKm = false;
            // convert km to m
            if constexpr(O3D) {
                for(int32_t iy = 0; iy < bdinfotb->NPts.y; ++iy) {
                    for(int32_t ix = 0; ix < bdinfotb->NPts.x; ++ix) {
                        vec3 &x = bdinfotb->bd[ix * bdinfotb->NPts.y + iy].x;
                        x.x *= RL(1000.0);
                        x.y *= RL(1000.0);
                    }
                }
            } else {
                for(int32_t i = 0; i < bdinfotb->NPts; ++i)
                    bdinfotb->bd[i].x[0] *= RL(1000.0);
            }
        }

        if(!bdinfotb->dirty) return;
        bdinfotb->dirty = false;

        ComputeBdryTangentNormal(params, bdinfotb);

        if constexpr(!O3D) {
            // convert range-dependent geoacoustic parameters from user to program units
            if(bdinfotb->type[1] == 'L') {
                for(int32_t iSeg = 0; iSeg < bdinfotb->NPts; ++iSeg) {
                    // compressional wave speed
                    bdinfotb->bd[iSeg].hs.cP = crci(
                        params, RL(1.0e20), bdinfotb->bd[iSeg].hs.alphaR,
                        bdinfotb->bd[iSeg].hs.alphaI, {'W', ' '});
                    // shear         wave speed
                    bdinfotb->bd[iSeg].hs.cS = crci(
                        params, RL(1.0e20), bdinfotb->bd[iSeg].hs.betaR,
                        bdinfotb->bd[iSeg].hs.betaI, {'W', ' '});
                }
            }
        }
    }

    virtual void Finalize(bhcParams<O3D> &params) const override
    {
        BdryInfoTopBot<O3D> *bdinfotb = GetBdryInfoTopBot(params);
        trackdeallocate(params, bdinfotb->bd);
    }

private:
    constexpr static real NegTop = ISTOP ? RL(-1.0) : RL(1.0);

    char &GetModeFlag(bhcParams<O3D> &params) const
    {
        if constexpr(ISTOP)
            return params.Bdry->Top.hs.Opt[4];
        else
            return params.Bdry->Bot.hs.Opt[1];
    }
    bool IsFile(bhcParams<O3D> &params) const
    {
        char flag = GetModeFlag(params);
        return flag == '~' || flag == '*';
    }
    real BdryDepth(bhcParams<O3D> &params) const
    {
        if constexpr(ISTOP)
            return params.Bdry->Top.hs.Depth;
        else
            return params.Bdry->Bot.hs.Depth;
    }
    BdryInfoTopBot<O3D> *GetBdryInfoTopBot(bhcParams<O3D> &params) const
    {
        if constexpr(ISTOP)
            return &params.bdinfo->top;
        else
            return &params.bdinfo->bot;
    }

    constexpr static const char *s_atibty              = ISTOP ? "ati" : "bty";
    constexpr static const char *s_ATIBTY              = ISTOP ? "ATI" : "BTY";
    constexpr static const char *s_altimetrybathymetry = ISTOP ? "altimetry"
                                                               : "bathymetry";
    constexpr static const char *s_AltimetryBathymetry = ISTOP ? "Altimetry"
                                                               : "Bathymetry";
    constexpr static const char *s_topbottom           = ISTOP ? "top" : "bottom";
    constexpr static const char *s_risesdrops          = ISTOP ? "rises above highest"
                                                               : "drops below lowest";

    /**
     * Does some pre-processing on the boundary points to pre-compute segment
     * lengths  (.Len),
     * tangents (.t, .nodet),
     * normals  (.n, .noden), and
     * curvatures (.kappa)
     */
    inline void ComputeBdryTangentNormal(
        const bhcParams<O3D> &params, BdryInfoTopBot<O3D> *bd) const
    {
        typename TmplInt12<O3D>::type NPts = bd->NPts;
        vec3 tvec;

        if constexpr(O3D) {
            // normals on triangle faces
            for(int32_t ix = 0; ix < NPts.x - 1; ++ix) {
                for(int32_t iy = 0; iy < NPts.y - 1; ++iy) {
                    // coordinates of corner nodes, moving counter-clockwise around the
                    // rectangle
                    vec3 p1 = bd->bd[(ix)*NPts.y + iy].x;
                    vec3 p2 = bd->bd[(ix + 1) * NPts.y + iy].x;
                    vec3 p3 = bd->bd[(ix + 1) * NPts.y + iy + 1].x;
                    vec3 p4 = bd->bd[(ix)*NPts.y + iy + 1].x;

                    // edges for triangle 1
                    vec3 u = p2 - p1; // tangent along one edge
                    vec3 v = p3 - p1; // tangent along another edge

                    // normal vector is the cross-product of the edge tangents
                    vec3 n1 = glm::cross(u, v);
                    if constexpr(ISTOP) n1 = -n1;

                    bd->bd[ix * NPts.y + iy].n1 = n1 / glm::length(n1); // scale to make
                                                                        // it a unit
                                                                        // normal

                    // edges for triangle 2
                    u = p3 - p1; // tangent along one edge
                    v = p4 - p1; // tangent along another edge

                    // normal vector is the cross-product of the edge tangents
                    vec3 n2 = glm::cross(u, v);
                    if constexpr(ISTOP) n2 = -n2;

                    bd->bd[ix * NPts.y + iy].n2 = n2 / glm::length(n2); // scale to make
                                                                        // it a unit
                                                                        // normal
                }
            }

            // normals at nodes
            // use forward, centered, or backward difference formulas
            for(int32_t ix = 0; ix < NPts.x; ++ix) {
                for(int32_t iy = 0; iy < NPts.y; ++iy) {
                    real mx, my;
                    if(ix == 0) {
                        mx = (bd->bd[(ix + 1) * NPts.y + iy].x.z
                              - bd->bd[(ix)*NPts.y + iy].x.z)
                            / (bd->bd[(ix + 1) * NPts.y + iy].x.x
                               - bd->bd[(ix)*NPts.y + iy].x.x);
                    } else if(ix == NPts.x - 1) {
                        mx = (bd->bd[(ix)*NPts.y + iy].x.z
                              - bd->bd[(ix - 1) * NPts.y + iy].x.z)
                            / (bd->bd[(ix)*NPts.y + iy].x.x
                               - bd->bd[(ix - 1) * NPts.y + iy].x.x);
                    } else {
                        mx = (bd->bd[(ix + 1) * NPts.y + iy].x.z
                              - bd->bd[(ix - 1) * NPts.y + iy].x.z)
                            / (bd->bd[(ix + 1) * NPts.y + iy].x.x
                               - bd->bd[(ix - 1) * NPts.y + iy].x.x);
                    }

                    if(iy == 0) {
                        my = (bd->bd[(ix)*NPts.y + iy + 1].x.z
                              - bd->bd[(ix)*NPts.y + iy].x.z)
                            / (bd->bd[(ix)*NPts.y + iy + 1].x.y
                               - bd->bd[(ix)*NPts.y + iy].x.y);
                    } else if(iy == NPts.y - 1) {
                        my = (bd->bd[(ix)*NPts.y + iy].x.z
                              - bd->bd[(ix)*NPts.y + iy - 1].x.z)
                            / (bd->bd[(ix)*NPts.y + iy].x.y
                               - bd->bd[(ix)*NPts.y + iy - 1].x.y);
                    } else {
                        my = (bd->bd[(ix)*NPts.y + iy + 1].x.z
                              - bd->bd[(ix)*NPts.y + iy - 1].x.z)
                            / (bd->bd[(ix)*NPts.y + iy + 1].x.y
                               - bd->bd[(ix)*NPts.y + iy - 1].x.y);
                    }

                    vec3 n = vec3(-mx, -my, RL(1.0)); // this is a normal to the surface

                    if(ix < NPts.x - 1 && iy < NPts.y - 1) {
                        // xx term
                        // this is the angle at each node
                        bd->bd[(ix)*NPts.y + iy].phi_xx = STD::atan2(n.z, n.x);

                        // xy term
                        tvec = bd->bd[(ix + 1) * NPts.y + iy + 1].x
                            - bd->bd[(ix)*NPts.y + iy].x;
                        real Len = STD::sqrt(SQ(tvec.x) + SQ(tvec.y));
                        tvec /= Len;
                        // this is the angle at each node
                        bd->bd[(ix)*NPts.y + iy].phi_xy
                            = STD::atan2(n.z, n.x * tvec.x + n.y * tvec.y);

                        // yy term
                        // this is the angle at each node
                        bd->bd[(ix)*NPts.y + iy].phi_yy = STD::atan2(n.z, n.y);
                    }

                    bd->bd[(ix)*NPts.y + iy].Noden_unscaled = n;
                    bd->bd[(ix)*NPts.y + iy].Noden          = n / glm::length(n);
                }
            }

        } else {
            // compute tangent and outward-pointing normal to each bottom segment
            // tBdry[0][:] = xBdry[0][1:NPts-1] - xBdry[0][0:NPts-2]
            // tBdry[1][:] = xBdry[1][1:NPts-1] - xBdry[1][0:NPts-2]
            // above caused compiler problems
            // LP: C++ obviously does not have vector slicing, but you get the idea.

            for(int32_t ii = 0; ii < NPts - 1; ++ii) {
                bd->bd[ii].t  = bd->bd[ii + 1].x - bd->bd[ii].x;
                bd->bd[ii].Dx = bd->bd[ii].t[1] / bd->bd[ii].t[0]; // first derivative
                // printf("Dx, t %g %g %g\n", bd->bd[ii].Dx, bd->bd[ii].x,
                //     (FL(1.0) / (bd->bd[ii].x[1] / FL(500.0)));

                // normalize the tangent vector
                bd->bd[ii].Len = glm::length(bd->bd[ii].t);
                bd->bd[ii].t /= bd->bd[ii].Len;

                bd->bd[ii].n[0] = -NegTop * bd->bd[ii].t[1];
                bd->bd[ii].n[1] = NegTop * bd->bd[ii].t[0];
            }
        }

        if(bd->type[0] == 'C') {
            // curvilinear option

            if constexpr(O3D) {
                // compute derivative as centered difference between two nodes
                // compute curvatures in each segment

                // - sign below because the node normal = vec3(-mx, -my, RL(1.0))
                for(int32_t ix = 0; ix < NPts.x - 1; ++ix) {
                    for(int32_t iy = 0; iy < NPts.y - 1; ++iy) {
                        real Len;

                        // z_xx (difference in x of z_x)
                        bd->bd[ix * NPts.y + iy].z_xx
                            = -(bd->bd[(ix + 1) * NPts.y + iy].Noden_unscaled.x
                                - bd->bd[(ix)*NPts.y + iy].Noden_unscaled.x)
                            / (bd->bd[(ix + 1) * NPts.y + iy].x.x
                               - bd->bd[(ix)*NPts.y + iy].x.x);

                        tvec = bd->bd[(ix + 1) * NPts.y + iy].x
                            - bd->bd[(ix)*NPts.y + iy].x;
                        Len = STD::sqrt(SQ(tvec.x) + SQ(tvec.z));
                        // this is curvature = dphi/ds
                        bd->bd[ix * NPts.y + iy].kappa_xx
                            = (bd->bd[(ix + 1) * NPts.y + iy].phi_xx
                               - bd->bd[(ix)*NPts.y + iy].phi_xx)
                            / Len;

                        // z_xy (difference in y of z_x)
                        bd->bd[ix * NPts.y + iy].z_xy
                            = -(bd->bd[(ix)*NPts.y + iy + 1].Noden_unscaled.x
                                - bd->bd[(ix)*NPts.y + iy].Noden_unscaled.x)
                            / (bd->bd[(ix)*NPts.y + iy + 1].x.y
                               - bd->bd[(ix)*NPts.y + iy].x.y);

                        // LP: This is overwritten by "new" below.
                        /*
                        tvec = bd->bd[(ix+1)*NPts.y+iy+1].x - bd->bd[(ix  )*NPts.y+iy ].x;
                        Len = glm::length(tvec);
                        // this is curvature = dphi/ds
                        bd->bd[ix*NPts.y+iy].kappa_xy = (bd->bd[(ix+1)*NPts.y+iy+1].phi_xy
                        - bd->bd[(ix  )*NPts.y+iy  ].phi_xy) / Len;
                        */

                        // new
                        tvec = bd->bd[(ix)*NPts.y + iy + 1].x
                            - bd->bd[(ix)*NPts.y + iy].x;
                        Len = STD::sqrt(SQ(tvec.y) + SQ(tvec.z));
                        // this is curvature = dphi/ds
                        bd->bd[ix * NPts.y + iy].kappa_xy
                            = (bd->bd[(ix)*NPts.y + iy + 1].phi_xx
                               - bd->bd[(ix)*NPts.y + iy].phi_xx)
                            / Len;

                        // z_yy (difference in y of z_y)
                        bd->bd[ix * NPts.y + iy].z_yy
                            = -(bd->bd[(ix)*NPts.y + iy + 1].Noden_unscaled.y
                                - bd->bd[(ix)*NPts.y + iy].Noden_unscaled.y)
                            / (bd->bd[(ix)*NPts.y + iy + 1].x.y
                               - bd->bd[(ix)*NPts.y + iy].x.y);

                        tvec = bd->bd[(ix)*NPts.y + iy + 1].x
                            - bd->bd[(ix)*NPts.y + iy].x;
                        Len = STD::sqrt(SQ(tvec.y) + SQ(tvec.z));
                        // this is curvature = dphi/ds
                        bd->bd[ix * NPts.y + iy].kappa_yy
                            = (bd->bd[(ix)*NPts.y + iy + 1].phi_xx
                               - bd->bd[(ix)*NPts.y + iy].phi_xx)
                            / Len;

                        // introduce Len factor per Eq. 4.4.18 in Cerveny's book
                        Len = glm::length(bd->bd[ix * NPts.y + iy].Noden_unscaled);
                        bd->bd[ix * NPts.y + iy].z_xx /= Len;
                        bd->bd[ix * NPts.y + iy].z_xy /= Len;
                        bd->bd[ix * NPts.y + iy].z_yy /= Len;
                    }
                }

                // LP: Last row and column data is uninitialized; make sure it is never
                // used.
                for(int32_t ix = 0; ix < NPts.x; ++ix) {
                    bd->bd[ix * NPts.y + (NPts.y - 1)].z_xx     = NAN;
                    bd->bd[ix * NPts.y + (NPts.y - 1)].z_xy     = NAN;
                    bd->bd[ix * NPts.y + (NPts.y - 1)].z_yy     = NAN;
                    bd->bd[ix * NPts.y + (NPts.y - 1)].kappa_xx = NAN;
                    bd->bd[ix * NPts.y + (NPts.y - 1)].kappa_xy = NAN;
                    bd->bd[ix * NPts.y + (NPts.y - 1)].kappa_yy = NAN;
                }
                for(int32_t iy = 0; iy < NPts.y; ++iy) {
                    bd->bd[(NPts.x - 1) * NPts.y + iy].z_xx     = NAN;
                    bd->bd[(NPts.x - 1) * NPts.y + iy].z_xy     = NAN;
                    bd->bd[(NPts.x - 1) * NPts.y + iy].z_yy     = NAN;
                    bd->bd[(NPts.x - 1) * NPts.y + iy].kappa_xx = NAN;
                    bd->bd[(NPts.x - 1) * NPts.y + iy].kappa_xy = NAN;
                    bd->bd[(NPts.x - 1) * NPts.y + iy].kappa_yy = NAN;
                }

            } else {
                // compute tangent and normal at node by averaging normals on adjacent
                // segments averaging two centered differences is equivalent to forming a
                // single centered difference of two steps ...
                for(int32_t ii = 1; ii < NPts - 1; ++ii) {
                    real sss = bd->bd[ii - 1].Len / (bd->bd[ii - 1].Len + bd->bd[ii].Len);
                    sss      = FL(0.5); // LP: BUG? Line above is overwritten.
                    bd->bd[ii].Nodet = (FL(1.0) - sss) * bd->bd[ii - 1].t
                        + sss * bd->bd[ii].t;
                }

                bd->bd[0].Nodet        = vec2(FL(1.0), FL(0.0)); // tangent left-end  node
                bd->bd[NPts - 1].Nodet = vec2(FL(1.0), FL(0.0)); // tangent right-end node

                for(int32_t ii = 0; ii < NPts; ++ii) {
                    bd->bd[ii].Noden[0] = -NegTop * bd->bd[ii].Nodet[1];
                    bd->bd[ii].Noden[1] = NegTop * bd->bd[ii].Nodet[0];
                }

                // compute curvature in each segment
                if(NPts < 2) {
                    EXTERR("Invalid value of NPts in ComputeBdryTangentNormal");
                }
                real curphi = STD::atan2(bd->bd[0].Nodet[1], bd->bd[0].Nodet[0]);
                for(int32_t ii = 0; ii < NPts - 1; ++ii) {
                    real nextphi
                        = STD::atan2(bd->bd[ii + 1].Nodet[1], bd->bd[ii + 1].Nodet[0]);
                    // this is curvature = dphi/ds
                    bd->bd[ii].kappa = (nextphi - curphi) / bd->bd[ii].Len;
                    // second derivative
                    bd->bd[ii].Dxx = (bd->bd[ii + 1].Dx - bd->bd[ii].Dx)
                        / (bd->bd[ii + 1].x[0] - bd->bd[ii].x[0]);
                    // derivative in direction of tangent
                    bd->bd[ii].Dss = bd->bd[ii].Dxx * CUBE(bd->bd[ii].t[0]);
                    // printf("kappa, Dss, Dxx %g %g %g %g %g %g %g\n", bd->bd[ii].kappa,
                    // bd->bd[ii].Dss, bd->bd[ii].Dxx,
                    //    FL(1.0) / ((FL(8.0) / SQ(FL(1000.0))) *
                    //    CUBE(STD::abs(bd->bd[ii].x[1]))), bd->bd[ii].x[1], FL(-1.0) /
                    //    (FL(4.0) * CUBE(bd->bd[ii].x[1]) / FL(1000000.0)),
                    //    bd->bd[ii].x[1]);

                    bd->bd[ii].kappa = bd->bd[ii].Dss; // over-ride kappa !!!!!
                    curphi           = nextphi;
                }
            }

        } else {
            if constexpr(O3D) {
                for(int32_t ix = 0; ix < NPts.x; ++ix) {
                    for(int32_t iy = 0; iy < NPts.y; ++iy) {
                        bd->bd[ix * NPts.y + iy].z_xx = RL(0.0);
                        bd->bd[ix * NPts.y + iy].z_xy = RL(0.0);
                        bd->bd[ix * NPts.y + iy].z_yy = RL(0.0);

                        bd->bd[ix * NPts.y + iy].kappa_xx = RL(0.0);
                        bd->bd[ix * NPts.y + iy].kappa_xy = RL(0.0);
                        bd->bd[ix * NPts.y + iy].kappa_yy = RL(0.0);
                    }
                }
            } else {
                for(int32_t i = 0; i < NPts; ++i) bd->bd[i].kappa = FL(0.0);
            }
        }
    }
};

template<bool O3D> using Altimetry  = Boundary<O3D, true>;
template<bool O3D> using Bathymetry = Boundary<O3D, false>;

}} // namespace bhc::module
