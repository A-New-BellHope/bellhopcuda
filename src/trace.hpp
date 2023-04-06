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
#include "common_run.hpp"
#include "step.hpp"
#include "reflect.hpp"
#include "boundary.hpp"
#include "influence.hpp"

namespace bhc {

/**
 * Calculates the distances to the boundaries
 * Formula differs from JKPS because code uses outward pointing normals
 *
 * rayx: ray coordinate
 * Topx, Botx: top, bottom coordinate
 * Topn, Botn: top, bottom normal vector (outward)
 * DistTop, DistBot: distance (normal to bdry) from the ray to top, bottom boundary
 */
template<bool O3D> HOST_DEVICE inline void Distances(
    const VEC23<O3D> &rayx, const VEC23<O3D> &Topx, const VEC23<O3D> &Botx,
    const VEC23<O3D> &Topn, const VEC23<O3D> &Botn, real &DistTop, real &DistBot)
{
    VEC23<O3D> dTop = rayx - Topx; // vector pointing from top    bdry to ray
    VEC23<O3D> dBot = rayx - Botx; // vector pointing from bottom bdry to ray
    DistTop         = -glm::dot(Topn, dTop);
    DistBot         = -glm::dot(Botn, dBot);
}

/**
 * LP: Not a typo that this is templated on O3D only. In BELLHOP3D, this happens
 * before the Nx2D/3D split, so it is only dependent on the ocean's dimensionality.
 */
template<typename CFG, bool O3D> HOST_DEVICE inline SSPOutputs<O3D> RayStartNominalSSP(
    int32_t isx, int32_t isy, int32_t isz, [[maybe_unused]] real alpha, SSPSegState &iSeg,
    const Position *Pos, const SSPStructure *ssp, ErrState *errState, VEC23<O3D> &xs,
    VEC23<O3D> &tinit)
{
    if constexpr(O3D) {
        xs    = vec3(Pos->Sx[isx], Pos->Sy[isy], Pos->Sz[isz]);
        tinit = vec3(FL(0.0), FL(0.0), FL(1.0));
    } else {
        xs    = vec2(FL(0.0), Pos->Sz[isz]); // x-y [LP: r-z] coordinate of the source
        tinit = vec2(STD::cos(alpha), STD::sin(alpha));
    }
    SSPOutputs<O3D> o;
    // LP: Not a typo that this is templated on O3D only; see function comment.
    EvaluateSSP<CFG, O3D, O3D>(xs, tinit, o, Origin<O3D, O3D>(), ssp, iSeg, errState);
    return o;
}

/**
 * LP: Pulled out ray update loop initialization. Returns whether to continue
 * with the ray trace. Only call for valid ialpha w.r.t. Angles->iSingleAlpha.
 * Original comments follow.
 *
 * DistBegTop etc.: Distances from ray beginning, end to top and bottom
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline bool RayInit(
    RayInitInfo &rinit, VEC23<O3D> &xs, rayPt<R3D> &point0, VEC23<O3D> &gradc,
    real &DistBegTop, real &DistBegBot, Origin<O3D, R3D> &org, SSPSegState &iSeg,
    BdryState<O3D> &bds, BdryType &Bdry, const BdryType *ConstBdry,
    const BdryInfo<O3D> *bdinfo, const SSPStructure *ssp, const Position *Pos,
    const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const SBPInfo *sbp, ErrState *errState)
{
    if(rinit.isz < 0 || rinit.isz >= Pos->NSz || rinit.ialpha < 0
       || rinit.ialpha >= Angles->alpha.n
       || (O3D
           && (rinit.isx < 0 || rinit.isx >= Pos->NSx || rinit.isy < 0
               || rinit.isy >= Pos->NSy || rinit.ibeta < 0
               || rinit.ibeta >= Angles->beta.n))) {
        RunError(errState, BHC_ERR_RAYINIT);
        // printf("Invalid ray init indexes!\n");
        return false;
    }
    // printf("Tracing azimuthal %d declination %d beam\n", rinit.ibeta, rinit.ialpha);

    // LP: This part from BellhopCore

    rinit.alpha        = Angles->alpha.angles[rinit.ialpha]; // initial angle
    rinit.SrcDeclAngle = RadDeg * rinit.alpha; // take-off declination angle in degrees
    if constexpr(O3D) {
        rinit.beta         = Angles->beta.angles[rinit.ibeta];
        rinit.SrcAzimAngle = RadDeg * rinit.beta; // take-off azimuthal   angle in degrees
    } else {
        rinit.beta = rinit.SrcAzimAngle = NAN;
    }

    iSeg.x = iSeg.y = iSeg.z = iSeg.r = 0;
    VEC23<O3D> tinit;
    SSPOutputs<O3D> o = RayStartNominalSSP<CFG, O3D>(
        rinit.isx, rinit.isy, rinit.isz, rinit.alpha, iSeg, Pos, ssp, errState, xs,
        tinit);
    gradc = o.gradc;

    if constexpr(O3D && !R3D) {
        org.xs      = xs;
        org.tradial = vec2(STD::cos(rinit.beta), STD::sin(rinit.beta));
    }

    if constexpr(!O3D) {
        // Are there enough beams?
        real DalphaOpt = STD::sqrt(
            o.ccpx.real() / (FL(6.0) * freqinfo->freq0 * Pos->Rr[Pos->NRr - 1]));
        int32_t NalphaOpt = 2
            + (int)((Angles->alpha.angles[Angles->alpha.n - 1] - Angles->alpha.angles[0])
                    / DalphaOpt);

        if(IsCoherentRun(Beam) && Angles->alpha.n < NalphaOpt && rinit.ialpha == 0) {
            RunWarning(errState, BHC_WARN_TOO_FEW_BEAMS);
            // printf(
            //     "Warning in " BHC_PROGRAMNAME
            //     " : Too few beams\nNalpha should be at least = %d\n",
            //     NalphaOpt);
        }
    }

    int32_t ibp = BinarySearchLEQ(sbp->SrcBmPat, sbp->NSBPPts, 2, 0, rinit.SrcDeclAngle);
    // LP: Our function won't ever go outside the table, but we need to limit it
    // to 2 from the end.
    ibp = bhc::min(ibp, sbp->NSBPPts - 2);

    // linear interpolation to get amplitude
    real s = (rinit.SrcDeclAngle - sbp->SrcBmPat[2 * ibp + 0])
        / (sbp->SrcBmPat[2 * (ibp + 1) + 0] - sbp->SrcBmPat[2 * ibp + 0]);
    float Amp0 = (FL(1.0) - s) * sbp->SrcBmPat[2 * ibp + 1]
        + s * sbp->SrcBmPat[2 * (ibp + 1) + 1]; // initial amplitude

    // Lloyd mirror pattern for semi-coherent option
    if(IsSemiCoherentRun(Beam)) {
        float omega = FL(2.0) * REAL_PI * freqinfo->freq0;
        Amp0 *= STD::sqrt(FL(2.0))
            * STD::abs(STD::sin(omega / o.ccpx.real() * DEP(xs) * STD::sin(rinit.alpha)));
    }

    // LP: This part from TraceRay

    VEC23<R3D> tinit2;
    if constexpr(R3D) {
        tinit2 = vec3(
            STD::cos(rinit.alpha) * STD::cos(rinit.beta),
            STD::cos(rinit.alpha) * STD::sin(rinit.beta), STD::sin(rinit.alpha));
    } else if constexpr(O3D) {
        tinit2 = vec2(STD::cos(rinit.alpha), STD::sin(rinit.alpha));
    } else {
        tinit2 = tinit;
    }
    if constexpr(O3D && !R3D) {
        point0.x = vec2(RL(0.0), xs.z);
    } else {
        point0.x = xs;
    }
    point0.c         = o.ccpx.real();
    point0.t         = tinit2 / o.ccpx.real();
    point0.tau       = cpx(FL(0.0), FL(0.0));
    point0.Amp       = Amp0;
    point0.Phase     = FL(0.0);
    point0.NumTopBnc = 0;
    point0.NumBotBnc = 0;
    if constexpr(R3D) {
        point0.p = mat2x2(FL(1.0)); // LP: identity
        point0.q = mat2x2(FL(0.0)); // LP: zero matrix
        // point0.DetQ  = epsilon1 * epsilon2 // LP: commented out
        point0.phi = FL(0.0);
    } else {
        point0.p = vec2(FL(1.0), FL(0.0));
        point0.q = vec2(FL(0.0), FL(1.0));
    }

    if constexpr(!O3D) {
        // second component of qv is not used in geometric beam tracing
        // set I.C. to 0 in hopes of saving run time
        // LP: I don't think modern CPUs / GPUs have early return from the
        // floating-point unit if one of the operands is zero, as the floating-
        // point unit is pipelined and produces results every clock anyway.
        // LP: Also, this only applies to geometric hat beams in Cartesian
        // coordinates specified as 'G' rather than '^', not any of the other
        // geometric types.
        if(Beam->RunType[1] == 'G') point0.q = vec2(FL(0.0), FL(0.0));
    }

    Bdry = *ConstBdry;
    if constexpr(O3D) {
        bds.top.Iseg.x = bds.top.Iseg.y = 0;
        bds.bot.Iseg.x = bds.bot.Iseg.y = 0;
    } else {
        bds.top.Iseg = bds.bot.Iseg = 0;
    }
    GetBdrySeg<O3D>(
        xs, RayToOceanT(point0.t, org), bds.top, &bdinfo->top, Bdry.Top, true, true,
        errState); // identify the top    segment above the source
    GetBdrySeg<O3D>(
        xs, RayToOceanT(point0.t, org), bds.bot, &bdinfo->bot, Bdry.Bot, false, true,
        errState); // identify the bottom segment below the source

    Distances<O3D>(
        xs, bds.top.x, bds.bot.x, bds.top.n, bds.bot.n, DistBegTop, DistBegBot);

    if(DistBegTop <= FL(0.0) || DistBegBot <= FL(0.0)) {
        RunWarning(errState, BHC_WARN_SOURCE_OUTSIDE_BOUNDARIES);
#ifdef STEP_DEBUGGING
        printf("Terminating the ray trace because the source is on or outside the "
               "boundaries\n");
        if(DistBegTop <= FL(0.0)) {
            if constexpr(R3D) {
                printf(
                    "point0.x %g,%g,%g bds.top.x %g,%g,%g bds.top.n %g,%g,%g DistBegTop "
                    "%g\n",
                    point0.x.x, point0.x.y, point0.x.z, bds.top.x.x, bds.top.x.y,
                    bds.top.x.z, bds.top.n.x, bds.top.n.y, bds.top.n.z, DistBegTop);
            } else if constexpr(O3D) {
                printf(
                    "point0.x %g,%g bds.top.x %g,%g,%g bds.top.n %g,%g,%g DistBegTop "
                    "%g\n",
                    point0.x.x, point0.x.y, bds.top.x.x, bds.top.x.y, bds.top.x.z,
                    bds.top.n.x, bds.top.n.y, bds.top.n.z, DistBegTop);
            } else {
                printf(
                    "point0.x %g,%g bds.top.x %g,%g bds.top.n %g,%g DistBegTop %g\n",
                    point0.x.x, point0.x.y, bds.top.x.x, bds.top.x.y, bds.top.n.x,
                    bds.top.n.y, DistBegTop);
            }
        } else {
            if constexpr(R3D) {
                printf(
                    "point0.x %g,%g,%g bds.bot.x %g,%g,%g bds.bot.n %g,%g,%g DistBegBot "
                    "%g\n",
                    point0.x.x, point0.x.y, point0.x.z, bds.bot.x.x, bds.bot.x.y,
                    bds.bot.x.z, bds.bot.n.x, bds.bot.n.y, bds.bot.n.z, DistBegBot);
            } else if constexpr(O3D) {
                printf(
                    "point0.x %g,%g bds.bot.x %g,%g,%g bds.bot.n %g,%g,%g DistBegBot "
                    "%g\n",
                    point0.x.x, point0.x.y, bds.bot.x.x, bds.bot.x.y, bds.bot.x.z,
                    bds.bot.n.x, bds.bot.n.y, bds.bot.n.z, DistBegBot);
            } else {
                printf(
                    "point0.x %g,%g bds.bot.x %g,%g bds.bot.n %g,%g DistBegBot %g\n",
                    point0.x.x, point0.x.y, bds.bot.x.x, bds.bot.x.y, bds.bot.n.x,
                    bds.bot.n.y, DistBegBot);
            }
        }
#endif
        return false; // source must be within the medium
    }

    return true;
}

/**
 * LP: Pulled out contents of ray update loop. Returns the number of ray steps
 * taken (i.e. normally 1, or 2 if reflected).
 *
 * Returns whether reflection happened and therefore a second step was taken.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline bool RayUpdate(
    const rayPt<R3D> &point0, rayPt<R3D> &point1, rayPt<R3D> &point2, real &DistEndTop,
    real &DistEndBot, int32_t &iSmallStepCtr, const Origin<O3D, R3D> &org,
    SSPSegState &iSeg, BdryState<O3D> &bds, BdryType &Bdry, const BdryInfo<O3D> *bdinfo,
    const ReflectionInfo *refl, const SSPStructure *ssp, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const VEC23<O3D> &xs, ErrState *errState)
{
    bool topRefl, botRefl;
    Step<CFG, O3D, R3D>(
        point0, point1, bds, Beam, xs, org, ssp, iSeg, errState, iSmallStepCtr, topRefl,
        botRefl);
    /*
    if(point0.x == point1.x){
        printf("Ray did not move from (%g,%g), bailing\n", point0.x.x, point0.x.y);
        bail();
    }
    */

    VEC23<O3D> x_o = RayToOceanX(point1.x, org);
    VEC23<O3D> t_o = RayToOceanT(point1.t, org);
    GetBdrySeg<O3D>(x_o, t_o, bds.top, &bdinfo->top, Bdry.Top, true, false, errState);
    GetBdrySeg<O3D>(x_o, t_o, bds.bot, &bdinfo->bot, Bdry.Bot, false, false, errState);

    // Reflections?
    // Tests that ray at step is is inside, and ray at step is+1 is outside
    // to detect only a crossing from inside to outside
    // DistBeg is the distance at point0, which is saved
    // DistEnd is the distance at point1, which needs to be calculated
    Distances<O3D>(
        x_o, bds.top.x, bds.bot.x, bds.top.n, bds.bot.n, DistEndTop, DistEndBot);

    // LP: Merging these cases is important for GPU performance.
    if(topRefl || botRefl) {
#ifdef STEP_DEBUGGING
        printf(topRefl ? "Top reflecting\n" : "Bottom reflecting\n");
#endif
        const BdryInfoTopBot<O3D> &bdi     = topRefl ? bdinfo->top : bdinfo->bot;
        const BdryStateTopBot<O3D> &bdstb  = topRefl ? bds.top : bds.bot;
        const HSInfo &hs                   = topRefl ? Bdry.Top.hs : Bdry.Bot.hs;
        const ReflectionInfoTopBot &refltb = topRefl ? refl->top : refl->bot;
        ReflCurvature<O3D> rcurv;
        VEC23<R3D> tInt(RL(0.0));
        VEC23<O3D> nInt;

        if constexpr(O3D) {
            // LP: FORTRAN actually checks if the whole string is just "C", not just the
            // first char
            if(bdi.type[0] == 'C') {
                real s1 = (x_o.x - bdstb.x.x) / (bdstb.lSeg.x.max - bdstb.lSeg.x.min);
                real s2 = (x_o.y - bdstb.x.y) / (bdstb.lSeg.y.max - bdstb.lSeg.y.min);
                real m1 = FL(1.0) - s1;
                real m2 = FL(1.0) - s2;

                BdryPtFull<true> *bd00
                    = &bdi.bd[(bdstb.Iseg.x) * bdi.NPts.y + bdstb.Iseg.y];
                BdryPtFull<true> *bd01
                    = &bdi.bd[(bdstb.Iseg.x) * bdi.NPts.y + bdstb.Iseg.y + 1];
                BdryPtFull<true> *bd10
                    = &bdi.bd[(bdstb.Iseg.x + 1) * bdi.NPts.y + bdstb.Iseg.y];
                BdryPtFull<true> *bd11
                    = &bdi.bd[(bdstb.Iseg.x + 1) * bdi.NPts.y + bdstb.Iseg.y + 1];

                nInt = bd00->Noden * m1 * m2 + bd10->Noden * s1 * m2
                    + bd11->Noden * s1 * s2 + bd01->Noden * m1 * s2;
                rcurv.z_xx = bd00->z_xx;
                rcurv.z_xy = bd00->z_xy;
                rcurv.z_yy = bd00->z_yy;

                rcurv.kappa_xx = bd00->kappa_xx;
                rcurv.kappa_xy = bd00->kappa_xy;
                rcurv.kappa_yy = bd00->kappa_yy;
            } else {
                nInt       = bdstb.n;
                rcurv.z_xx = rcurv.z_xy = rcurv.z_yy = FL(0.0);
                rcurv.kappa_xx = rcurv.kappa_xy = rcurv.kappa_yy = FL(0.0);
            }
        } else {
            BdryPtFull<false> *bd0 = &bdi.bd[bdstb.Iseg];
            BdryPtFull<false> *bd1 = &bd0[1]; // LP: next segment
            // LP: FORTRAN actually checks if the whole string is just "C", not just the
            // first char
            if(bdi.type[0] == 'C') {
                real sss = glm::dot(point1.x - bdstb.x, bd0->t)
                    / bd0->Len; // proportional
                                // distance
                                // along
                                // segment
                nInt = (FL(1.0) - sss) * bd0->Noden + sss * bd1->Noden;
                tInt = (FL(1.0) - sss) * bd0->Nodet + sss * bd1->Nodet;
            } else {
                nInt = bd0->n; // normal is constant in a segment
                tInt = bd0->t;
            }
            rcurv.kappa = bd0->kappa;
        }

        Reflect<CFG, O3D, R3D>(
            point1, point2, hs, topRefl, tInt, nInt, rcurv, freqinfo->freq0, refltb, Beam,
            org, ssp, iSeg, errState);
        // Incrementing bounce count moved to Reflect
        x_o = RayToOceanX(point2.x, org);
        Distances<O3D>(
            x_o, bds.top.x, bds.bot.x, bds.top.n, bds.bot.n, DistEndTop, DistEndBot);
        return true;
    }

    return false;
}

/**
 * Has the ray left the box, lost its energy, escaped the boundaries, or
 * exceeded storage limit?
 * [LP: Nx2D only:]
 * this should be modified to have a single box
 * no need to test point.x.x, for instance, against several limits; calculate one limit in
 * advance LP: Also updates DistBegTop, DistBegBot.
 */
template<bool O3D, bool R3D> HOST_DEVICE inline bool RayTerminate(
    const rayPt<R3D> &point, int32_t &Nsteps, int32_t is, const VEC23<O3D> &xs,
    const int32_t &iSmallStepCtr, real &DistBegTop, real &DistBegBot,
    const real &DistEndTop, const real &DistEndBot, int32_t MaxPointsPerRay,
    const Origin<O3D, R3D> &org, [[maybe_unused]] const BdryInfo<O3D> *bdinfo,
    const BeamStructure<O3D> *Beam, ErrState *errState)
{
    bool leftbox, escapedboundaries, toomanysmallsteps;
    if constexpr(O3D) {
        vec3 x_o = RayToOceanX(point.x, org);
        vec3 t_o = RayToOceanT(point.t, org);
        // LP: This condition was inexplicably commented out in 2022 revision of Nx2D,
        // see Fortran version readme.
        leftbox = IsOutsideBeamBoxDim<true, 0>(x_o, Beam, xs)
            || IsOutsideBeamBoxDim<true, 1>(x_o, Beam, xs)
            || IsOutsideBeamBoxDim<true, 2>(x_o, Beam, xs);
        real minx = bhc::max(bdinfo->bot.bd[0].x.x, bdinfo->top.bd[0].x.x);
        real miny = bhc::max(bdinfo->bot.bd[0].x.y, bdinfo->top.bd[0].x.y);
        real maxx = bhc::min(
            bdinfo->bot.bd[(bdinfo->bot.NPts.x - 1) * bdinfo->bot.NPts.y].x.x,
            bdinfo->top.bd[(bdinfo->top.NPts.x - 1) * bdinfo->top.NPts.y].x.x);
        real maxy = bhc::min(
            bdinfo->bot.bd[bdinfo->bot.NPts.y - 1].x.y,
            bdinfo->top.bd[bdinfo->top.NPts.y - 1].x.y);
        bool escaped0bdry, escapedNbdry;
        escaped0bdry = x_o.x < minx || x_o.y < miny;
        escapedNbdry = x_o.x > maxx || x_o.y > maxy;
        // LP: See discussion in GetBdrySeg for why this was changed from the
        // original BELLHOP3D version.
        escaped0bdry = escaped0bdry || (x_o.x == minx && t_o.x < RL(0.0))
            || (x_o.y == miny && t_o.y < RL(0.0));
        escapedNbdry = escapedNbdry || (x_o.x == maxx && t_o.x > RL(0.0))
            || (x_o.y == maxy && t_o.y > RL(0.0));
        escapedboundaries = escaped0bdry || escapedNbdry;
        toomanysmallsteps = iSmallStepCtr > 50;
    } else {
        leftbox = IsOutsideBeamBoxDim<false, 0>(point.x, Beam, xs)
            || IsOutsideBeamBoxDim<false, 1>(point.x, Beam, xs);
        escapedboundaries = (DistBegTop < FL(0.0) && DistEndTop < FL(0.0))
            || (DistBegBot < FL(0.0) && DistEndBot < FL(0.0));
        toomanysmallsteps = false; // LP: The small step counter is never checked in 2D.
    }
    bool lostenergy = point.Amp < FL(0.005);
    bool backward   = false;
    if constexpr(O3D && !R3D) {
        // backward = point.t.x < FL(0.0); // kills off a backward traveling ray
        // LP: Condition above is now (2022) commented out in Nx2D as well.
    }
    if(leftbox || lostenergy || escapedboundaries || backward || toomanysmallsteps) {
#ifdef STEP_DEBUGGING
        if(leftbox) {
            if constexpr(O3D) {
                printf(
                    "Ray left beam box (%g,%g,%g)\n", Beam->Box.x, Beam->Box.y,
                    Beam->Box.z);
            } else {
                printf("Ray left beam box (%g,%g)\n", Beam->Box.x, Beam->Box.y);
            }
        } else if(escapedboundaries) {
            printf(
                "Ray escaped boundaries DistBegTop %g DistEndTop %g DistBegBot %g "
                "DistEndBot %g\n",
                DistBegTop, DistEndTop, DistBegBot, DistEndBot);
        } else if(lostenergy) {
            printf("Ray energy dropped to %g\n", point.Amp);
        } else if(backward) {
            printf("Ray is going backwards\n");
        } else if(toomanysmallsteps) {
            printf("Too many small steps\n");
        } else {
            printf("Internal error in RayTerminate\n");
            bail();
        }
#endif
        Nsteps = is + 1;
        return true;
    } else if(is >= MaxPointsPerRay - 3) {
        RunWarning(errState, BHC_WARN_ONERAY_OUTOFMEMORY);
        // printf("Warning in TraceRay: Insufficient storage for ray trajectory\n");
        Nsteps = is;
        return true;
    }

    DistBegTop = DistEndTop;
    DistBegBot = DistEndBot;
    return false;
}

/**
 * Main ray tracing function for ray path output mode.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void MainRayMode(
    RayInitInfo &rinit, rayPt<R3D> *ray, int32_t &Nsteps, int32_t MaxPointsPerRay,
    Origin<O3D, R3D> &org, const BdryType *ConstBdry, const BdryInfo<O3D> *bdinfo,
    const ReflectionInfo *refl, const SSPStructure *ssp, const Position *Pos,
    const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const SBPInfo *sbp, ErrState *errState)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;

    if(!RayInit<CFG, O3D, R3D>(
           rinit, xs, ray[0], gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, sbp, errState)) {
        Nsteps = 1;
        return;
    }

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray

    while(true) {
        if(HasErrored(errState)) break;
        bool twoSteps = RayUpdate<CFG, O3D, R3D>(
            ray[is], ray[is + 1], ray[is + 2], DistEndTop, DistEndBot, iSmallStepCtr, org,
            iSeg, bds, Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs, errState);
        if(Nsteps >= 0 && is >= Nsteps) {
            Nsteps = is + 2;
            break;
        }
        is += (twoSteps ? 2 : 1);
        if(RayTerminate<O3D, R3D>(
               ray[is], Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, MaxPointsPerRay, org, bdinfo, Beam, errState))
            break;
    }
}

/**
 * Main ray tracing function for TL, eigen, and arrivals runs.
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void MainFieldModes(
    RayInitInfo &rinit, cpxf *uAllSources, const BdryType *ConstBdry,
    const BdryInfo<O3D> *bdinfo, const ReflectionInfo *refl, const SSPStructure *ssp,
    const Position *Pos, const AnglesStructure *Angles, const FreqInfo *freqinfo,
    const BeamStructure<O3D> *Beam, const SBPInfo *sbp, EigenInfo *eigen,
    const ArrInfo *arrinfo, ErrState *errState)
{
    real DistBegTop, DistEndTop, DistBegBot, DistEndBot;
    SSPSegState iSeg;
    VEC23<O3D> xs, gradc;
    BdryState<O3D> bds;
    BdryType Bdry;
    Origin<O3D, R3D> org;

    rayPt<R3D> point0, point1, point2;
    point2.c = NAN; // Silence incorrect g++ warning about maybe uninitialized;
    // it is always set when doing two steps, and not used otherwise
    InfluenceRayInfo<R3D> inflray;

    if(!RayInit<CFG, O3D, R3D>(
           rinit, xs, point0, gradc, DistBegTop, DistBegBot, org, iSeg, bds, Bdry,
           ConstBdry, bdinfo, ssp, Pos, Angles, freqinfo, Beam, sbp, errState)) {
        return;
    }

    Init_Influence<CFG, O3D, R3D>(
        inflray, point0, rinit, gradc, Pos, org, ssp, iSeg, Angles, freqinfo, Beam,
        errState);

    int32_t iSmallStepCtr = 0;
    int32_t is            = 0; // index for a step along the ray
    int32_t Nsteps        = 0; // not actually needed in TL mode, debugging only

    while(true) {
        if(HasErrored(errState)) break;
        bool twoSteps = RayUpdate<CFG, O3D, R3D>(
            point0, point1, point2, DistEndTop, DistEndBot, iSmallStepCtr, org, iSeg, bds,
            Bdry, bdinfo, refl, ssp, freqinfo, Beam, xs, errState);
        if(!Step_Influence<CFG, O3D, R3D>(
               point0, point1, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg, Pos,
               Beam, eigen, arrinfo, errState)) {
#ifdef STEP_DEBUGGING
            printf("Step_Influence terminated ray\n");
#endif
            break;
        }
        ++is;
        if(twoSteps) {
            if(!Step_Influence<CFG, O3D, R3D>(
                   point1, point2, inflray, is, uAllSources, ConstBdry, org, ssp, iSeg,
                   Pos, Beam, eigen, arrinfo, errState))
                break;
            point0 = point2;
            ++is;
        } else {
            point0 = point1;
        }
        if(RayTerminate<O3D, R3D>(
               point0, Nsteps, is, xs, iSmallStepCtr, DistBegTop, DistBegBot, DistEndTop,
               DistEndBot, MaxN, org, bdinfo, Beam, errState))
            break;
    }

    // printf("Nsteps %d\n", Nsteps);
}

} // namespace bhc
