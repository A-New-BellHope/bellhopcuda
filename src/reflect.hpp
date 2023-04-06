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

namespace bhc {

/**
 * Given an angle RInt%ThetaInt, returns the magnitude and
 * phase of the reflection coefficient (RInt%R, RInt%phi).
 *
 * Uses linear interpolation using the two nearest abscissas
 * Assumes phi has been unwrapped so that it varies smoothly.
 * I tried modifying it to allow a complex angle of incidence but
 * stopped when I realized I needed to fuss with a complex atan2 routine
 * LP: C++ supports complex atan, though there is no complex atan2.
 *
 * RInt: interpolated value of refl. coef.
 * rtb: Reflection coefficient table
 * NPts: # pts in refl. coef.
 */
HOST_DEVICE inline void InterpolateReflectionCoefficient(
    ReflectionCoef &RInt, const ReflectionInfoTopBot &rtb, ErrState *errState)
{
    int32_t iLeft, iRight, iMid;
    real alpha, thetaIntr;

    iLeft  = 0;
    iRight = rtb.NPts - 1;

    // LP: This was originally the FORTRAN version of RInt.theta.real(), but
    // theta is definitely already a real (originally double).
    thetaIntr = RInt.theta; // This should be unnecessary? probably used when I was doing
                            // complex angles

    // Three cases: ThetaInt left, in, or right of tabulated interval

    if(thetaIntr < rtb.r[iLeft].theta) {
        // iRight = 1;
        RInt.r   = FL(0.0); // rtb.r[iLeft].r
        RInt.phi = FL(0.0); // rtb.r[iLeft].phi
        RunWarning(errState, BHC_WARN_OUTSIDE_REFLCOEF);
        // printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
        //        "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
        //        thetaIntr, rtb.r[iLeft].theta);
    } else if(thetaIntr > rtb.r[iRight].theta) {
        // iLeft = rtb.NPts - 2;
        RInt.r   = FL(0.0); // rtb.r[iRight].r
        RInt.phi = FL(0.0); // rtb.r[iRight].phi
        // LP: The warning is commented out in BELLHOP in this case.
        // RunWarning(errState, BHC_WARN_OUTSIDE_REFLCOEF);
        // printf("Warning in InterpolateReflectionCoefficient : Refl. Coef. being "
        //        "set to 0 outside tabulated domain : angle = %f, lower limit = %f",
        //        thetaIntr, rtb.r[iRight].theta);
    } else {
        // Search for bracketing abscissas: STD::log2(rtb.NPts) stabs required for a
        // bracket

        while(iLeft != iRight - 1) {
            iMid = (iLeft + iRight) / 2;
            if(rtb.r[iMid].theta > thetaIntr) {
                iRight = iMid;
            } else {
                iLeft = iMid;
            }
        }

        // Linear interpolation for reflection coef

        alpha = (RInt.theta - rtb.r[iLeft].theta)
            / (rtb.r[iRight].theta - rtb.r[iLeft].theta);
        RInt.r   = (FL(1.0) - alpha) * rtb.r[iLeft].r + alpha * rtb.r[iRight].r;
        RInt.phi = (FL(1.0) - alpha) * rtb.r[iLeft].phi + alpha * rtb.r[iRight].phi;
    }
}

template<bool O3D, bool R3D> HOST_DEVICE inline ReflCurvature<R3D> OceanToRayCurvature(
    const ReflCurvature<O3D> &rcurv, const Origin<O3D, R3D> &org,
    [[maybe_unused]] bool isTop)
{
    static_assert(O3D || !R3D, "2D ocean but 3D rays not allowed!");
    if constexpr(O3D && !R3D) {
        ReflCurvature<false> rcurv_out;
        // mbp: use kappa_xx or z_xx?
        // rcurv_out.kappa = rcurv.z_xx * SQ(org.tradial.x)
        //     + FL(2.0) * rcurv.z_xy * org.tradial.x * org.tradial.y
        //     + rcurv.z_yy * SQ(org.tradial.y);
        rcurv_out.kappa = rcurv.kappa_xx * SQ(org.tradial.x)
            + FL(2.0) * rcurv.kappa_xy * org.tradial.x * org.tradial.y
            + rcurv.kappa_yy * SQ(org.tradial.y);
        if(isTop) rcurv_out.kappa = -rcurv_out.kappa;
        return rcurv_out;
    } else {
        return rcurv;
    }
}

/**
 * LP: Given that a reflection has occurred, reflect the ray/beam off the top
 * or bottom boundary.
 *
 * hs: half-space properties
 * isTop: Flag indicating bottom or top reflection
 * tBdry, nBdry: Tangent and normal to the boundary
 * rcurv: Boundary curvature
 * rtb: Reflection coefficient table
 */
template<typename CFG, bool O3D, bool R3D> HOST_DEVICE inline void Reflect(
    const rayPt<R3D> &oldPoint, rayPt<R3D> &newPoint, const HSInfo &hs, bool isTop,
    VEC23<R3D> tBdry, const VEC23<O3D> &nBdry, const ReflCurvature<O3D> &rcurv, real freq,
    const ReflectionInfoTopBot &rtb, [[maybe_unused]] const BeamStructure<O3D> *Beam,
    const Origin<O3D, R3D> &org, const SSPStructure *ssp, SSPSegState &iSeg,
    ErrState *errState)
{
    VEC23<R3D> nBdry_ray = OceanToRayT<O3D, R3D>(nBdry, org);
    if constexpr(O3D && !R3D) nBdry_ray *= RL(1.0) / glm::length(nBdry_ray);

    real Th = glm::dot(oldPoint.t, nBdry_ray); // component of ray tangent, normal to
                                               // boundary

    if constexpr(O3D) {
        // LP: tBdry is computed here for Nx2D and 3D. It is not only precomputed,
        // but computed differently from this formula, in 2D.
        tBdry = oldPoint.t - Th * nBdry_ray;
        tBdry *= RL(1.0) / glm::length(tBdry);
        // mbp: could also calculate tBdry as +/- of vec2(nBdry.y, -nBdry.x), but need
        // sign
    }

    real Tg = glm::dot(oldPoint.t, tBdry); // component of ray tangent, along boundary

    ReflCurvature<R3D> rcurv_ray = OceanToRayCurvature<O3D, R3D>(rcurv, org, isTop);

    // LP: Incrementing bounce count moved here
    newPoint.NumTopBnc = oldPoint.NumTopBnc + (int)isTop;
    newPoint.NumBotBnc = oldPoint.NumBotBnc + (1 - (int)isTop);
    newPoint.x         = oldPoint.x;
    newPoint.t = oldPoint.t - FL(2.0) * Th * nBdry_ray; // changing the ray direction

    // Calculate the change in curvature
    // Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

    // just to get c [LP: o.ccpx.real(); also, this comment is wrong, it is also using
    // o.gradc]
    SSPOutputs<R3D> o;
    EvaluateSSP<CFG, O3D, R3D>(
        oldPoint.x, (O3D ? newPoint.t : oldPoint.t), o, org, ssp, iSeg, errState);

    newPoint.c   = o.ccpx.real();
    newPoint.tau = oldPoint.tau;

    if constexpr(R3D) {
        vec3 rayt, rayn1, rayn2, rayt_tilde, rayn1_tilde, rayn2_tilde;
        CalcTangent_Normals(
            oldPoint, o.ccpx.real(), nBdry, rayt, rayn1, rayn2, RL(-1.0)); // incident
        CalcTangent_Normals(
            newPoint, o.ccpx.real(), nBdry, rayt_tilde, rayn1_tilde, rayn2_tilde,
            RL(-1.0)); // reflected

        // printf("point0 (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
        //     rayt.x, rayt.y, rayt.z, rayn1.x, rayn1.y, rayn1.z, rayn2.x, rayn2.y,
        //     rayn2.z);
        // printf("point1 (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
        //     rayt_tilde.x, rayt_tilde.y, rayt_tilde.z, rayn1_tilde.x, rayn1_tilde.y,
        //     rayn1_tilde.z, rayn2_tilde.x, rayn2_tilde.y, rayn2_tilde.z);

        // rotation matrix to get surface curvature in and perpendicular to the reflection
        // plane we use only the first two elements of the vectors because we want the
        // projection in the x-y plane
        vec2 t_rot = XYCOMP(rayt);
        t_rot *= RL(1.0) / glm::length(t_rot);
        vec2 n_rot = XYCOMP(rayn2);
        n_rot *= RL(1.0) / glm::length(n_rot);

        mat2x2 RotMat;
        RotMat[0] = t_rot; // LP: not a typo, the Fortran sets the columns
        RotMat[1] = n_rot;

        mat2x2 kappaMat; // LP: not clear why divided by 2 and then multiplied by 2
        kappaMat[0][0] = rcurv_ray.z_xx / FL(2.0);
        kappaMat[1][0] = rcurv_ray.z_xy / FL(2.0);
        kappaMat[0][1] = rcurv_ray.z_xy / FL(2.0);
        kappaMat[1][1] = rcurv_ray.z_yy / FL(2.0);

        // apply the rotation to get the matrix D of curvatures (see Popov 1977 for
        // definition of DMat) DMat = RotMat^T * kappaMat * RotMat, with RotMat
        // anti-symmetric
        mat2x2 DMatTemp = (FL(2.0) * kappaMat) * RotMat;
        mat2x2 DMat     = glm::transpose(RotMat) * DMatTemp;
        // PrintMatrix(DMat, "DMat");

        // normal and tangential derivatives of the sound speed
        real cn1jump = glm::dot(o.gradc, -rayn1_tilde - rayn1);
        real cn2jump = glm::dot(o.gradc, -rayn2_tilde - rayn2);
        real csjump  = -glm::dot(o.gradc, rayt_tilde - rayt);
        // printf("cn1jump cn2jump csjump %g, %g, %g\n", cn1jump, cn2jump, csjump);

        // // not sure if cn2 needs a sign flip also
        // if(isTop){
        //     cn1jump = -cn1jump; // flip sign for top reflection
        //     cn2jump = -cn2jump; // flip sign for top reflection
        // }

        vec3 e1, e2;
        RayNormal(oldPoint.t, oldPoint.phi, oldPoint.c, e1, e2);

        // LP: The routine below modifies newPoint in place, so have to copy first.
        newPoint.p   = oldPoint.p;
        newPoint.q   = oldPoint.q;
        newPoint.phi = oldPoint.phi;
        CurvatureCorrection3D<true>(
            newPoint, DMat, Tg, Th, cn1jump, cn2jump, csjump, rayn1, rayn2, e1, e2);
    } else {
        // incident unit ray tangent and normal
        VEC23<R3D> rayt = o.ccpx.real() * oldPoint.t; // unit tangent to ray
        VEC23<R3D> rayn = vec2(-rayt.y, rayt.x);      // unit normal  to ray

        // reflected unit ray tangent and normal (the reflected tangent, normal system has
        // a different orientation)
        VEC23<R3D> rayt_tilde = o.ccpx.real() * newPoint.t;         // unit tangent to ray
        VEC23<R3D> rayn_tilde = -vec2(-rayt_tilde.y, rayt_tilde.x); // unit normal  to ray
        // boundary curvature correction
        real rn = FL(2.0) * rcurv_ray.kappa / SQ(o.ccpx.real()) / Th;

        // get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th *
        // nbdry
        real cnjump = -glm::dot(o.gradc, rayn_tilde - rayn);
        real csjump = -glm::dot(o.gradc, rayt_tilde - rayt);

        if(isTop) {
            cnjump = -cnjump; // this is because the (t,n) system of the top boundary has
                              // a different sense to the bottom boundary
            rn = -rn;
        }

        real rm = Tg / Th; // this is tan( alpha ) where alpha is the angle of incidence
        rn += rm * (FL(2.0) * cnjump - rm * csjump) / SQ(o.ccpx.real());

        if(Beam->Type[2] == 'D') {
            rn *= FL(2.0);
        } else if(Beam->Type[2] == 'Z') {
            rn = FL(0.0);
        }

        newPoint.p = oldPoint.p + oldPoint.q * rn;
        newPoint.q = oldPoint.q;
    }

    // account for amplitude and phase change

    if(hs.bc == 'R') { // rigid
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase;
    } else if(hs.bc == 'V') { // vacuum
        newPoint.Amp   = oldPoint.Amp;
        newPoint.Phase = oldPoint.Phase + REAL_PI;
    } else if(hs.bc == 'F') { // file
        ReflectionCoef RInt;
        RInt.theta = RadDeg * STD::abs(STD::atan2(Th, Tg)); // angle of incidence
                                                            // (relative to normal to
                                                            // bathymetry)
        if(RInt.theta > FL(90.0))
            RInt.theta = FL(180.0) - RInt.theta; // reflection coefficient is symmetric
                                                 // about 90 degrees
        InterpolateReflectionCoefficient(RInt, rtb, errState);
        newPoint.Amp   = oldPoint.Amp * RInt.r;
        newPoint.Phase = oldPoint.Phase + RInt.phi;
    } else if(hs.bc == 'A' || hs.bc == 'G') { // half-space
        real omega = FL(2.0) * REAL_PI * freq;
        cpx Refl;
        if constexpr(O3D) {
            cpx gk = omega * Tg; // wavenumber in direction parallel to bathymetry
            // MINPOS prevents g95 [LP: gfortran] giving -zero, and wrong branch cut
            cpx gamma1Sq = SQ(omega / o.ccpx.real()) - SQ(gk) - J * REAL_MINPOS;
            cpx gamma2Sq = SQ(omega / hs.cP) - SQ(gk) - J * REAL_MINPOS;
            cpx gamma1   = STD::sqrt(-gamma1Sq);
            cpx gamma2   = STD::sqrt(-gamma2Sq);

            Refl = (hs.rho * gamma1 - o.rho * gamma2)
                / (hs.rho * gamma1 + o.rho * gamma2);

            // printf("%f %f %f %f %f\n", STD::abs(Refl), o.ccpx.real(), hs.cP, o.rho,
            // hs.rho);
            // LP: Removed this code as it was obviously modified to do nothing
            // (abs(anything) < 0 is always false), and users would likely be
            // unpleasantly surprised if it did do what it was intended to do.
            /*
            if constexpr(R3D) {
                // Hack to make a wall (where the bottom slope is more than 80 degrees) be
                // a perfect reflector
                if(STD::abs(
                       RadDeg
                       * STD::atan2(nBdry.z, glm::length(bhc::vec2(nBdry.x, nBdry.y))))
                   < 0) { // was 60 degrees
                    Refl = FL(1.0);
                }
            }
            */
        } else {
            cpx kx = omega * Tg; // wavenumber in direction parallel      to bathymetry
            cpx kz = omega * Th; // wavenumber in direction perpendicular to bathymetry
                                 // (in ocean)
            cpx kx2 = SQ(kx);

            // notation below is a bit misleading
            // kzS, kzP is really what I called gamma in other codes, and differs by a
            // factor of +/- i
            cpx f, g;
            if(hs.cS.real() > FL(0.0)) {
                cpx kzS2 = kx2 - SQ(omega / hs.cS);
                cpx kzP2 = kx2 - SQ(omega / hs.cP);
                cpx kzS  = STD::sqrt(kzS2);
                cpx kzP  = STD::sqrt(kzP2);
                cpx mu   = hs.rho * SQ(hs.cS);

                cpx y2 = (SQ(kzS2 + kx2) - RL(4.0) * kzS * kzP * kx2) * mu;
                cpx y4 = kzP * (kx2 - kzS2);

                f = SQ(omega) * y4;
                g = y2;
            } else {
                cpx kzP = STD::sqrt(kx2 - SQ(omega / hs.cP));

                // Intel and GFortran compilers return different branches of the SQRT for
                // negative reals LP: looks like this just means we want the positive
                // branch
                if(kzP.real() == RL(0.0) && kzP.imag() < RL(0.0)) kzP = -kzP;
                f = kzP;
                g = hs.rho;
            }

            Refl = -(o.rho * f - J * kz * g) / (o.rho * f + J * kz * g); // complex
                                                                         // reflection
                                                                         // coef.
            /*
            printf("cS cP rho (%g,%g) (%g,%g) %g\n", hs.cS.real(), hs.cS.imag(),
                hs.cP.real(), hs.cP.imag(), hs.rho);
            printf("kx kz f g Refl (%g,%g) (%g,%g) (%g,%g) (%g,%g) (%g,%g)\n",
                kx.real(), kx.imag(), kz.real(), kz.imag(), f.real(), f.imag(),
                g.real(), g.imag(), Refl.real(), Refl.imag());
            */
        }

        if(STD::abs(Refl) < FL(1.0e-5)) { // kill a ray that has lost its energy in
                                          // reflection
            newPoint.Amp   = FL(0.0);
            newPoint.Phase = oldPoint.Phase;
        } else {
            newPoint.Amp   = STD::abs(Refl) * oldPoint.Amp;
            newPoint.Phase = oldPoint.Phase + STD::atan2(Refl.imag(), Refl.real());

            if constexpr(!R3D) {
                // compute beam-displacement Tindle, Eq. (14)
                // needs a correction to beam-width as well ...
                // LP: most of these variables don't exist, likely very old code
                // if(kz2Sq.real() < FL(0.0)){
                //     rhoW = FL(1.0); // density of water
                //     rhoWSq = rhoW * rhoW;
                //     rhoHSSq = rhoHS * rhoHS;
                //     delta = FL(2.0) * gk * rhoW * rhoS * (kz1Sq - kz2Sq) /
                //         (kz1 * i * kz2 *
                //             (-rhoWSq * kz2Sq + rhoHSSq * kz1Sq));
                //     rv[is+1] = rv[is+1] + delta;
                // }

                if(Beam->Type[3] == 'S') { // beam displacement & width change (Seongil's
                                           // version)
                    cpx ch  = oldPoint.c / STD::conj(hs.cP);
                    real co = oldPoint.t.x * oldPoint.c;
                    real si = oldPoint.t.y * oldPoint.c;
                    real ck = omega / oldPoint.c;

                    cpx a    = FL(2.0) * hs.rho * (FL(1.0) - SQ(ch));
                    cpx b    = SQ(co) - SQ(ch);
                    cpx d    = SQ(hs.rho) * SQ(si) + b;
                    cpx sb   = STD::sqrt(b);
                    real cco = SQ(co);
                    real ssi = SQ(si);

                    cpx delta;
                    if(si != FL(0.0)) {
                        delta = a * co / si / (ck * sb * d); // Do we need an abs() on
                                                             // this???
                    } else {
                        delta = FL(0.0);
                    }

                    real pdelta = delta.real() / (oldPoint.c / co);
                    // LP: The spacing in the original version of this formula,
                    // the fact that several terms could be factored out to reduce
                    // computation, and the repeated divisons, lead me to believe
                    // that it may not be correct.
                    // Here is the original version with the weird spacing:
                    // ddelta = -a / (ck*sb*d) - a*cco / ssi / (ck*sb*d) + a*cco /
                    // (ck*b*sb*d)
                    //     -a*co / si / (ck*sb*d*d) * (FL(2.0)* SQ(hs.rho)
                    //     *si*co-FL(2.0)*co*si);
                    // Here is a version with things factored better:
                    cpx cksbd  = ck * sb * d;
                    cpx ddelta = a
                        * (cco / (cksbd * b) - (RL(1.0) + (cco / ssi)) / cksbd
                           - FL(2.0) * SQ(co) * (SQ(hs.rho) - RL(1.0)) / (cksbd * d));
                    real rddelta = -ddelta.real();
                    real sddelta = rddelta / STD::abs(rddelta);

                    if constexpr(!O3D) {
                        // next 3 lines have an update by Diana McCammon to allow a
                        // sloping bottom I think the formulas are good, but this won't be
                        // reliable because it doesn't have the logic that tracks crossing
                        // into new segments after the ray displacement.

                        real theta_bot = STD::atan(tBdry.y / tBdry.x); // bottom angle
                        newPoint.x.x   = newPoint.x.x
                            + delta.real() * STD::cos(theta_bot); // range displacement
                        newPoint.x.y = newPoint.x.y
                            + delta.real() * STD::sin(theta_bot); // depth displacement
                    } else {
                        newPoint.x.x = newPoint.x.x + delta.real(); // displacement
                    }
                    newPoint.tau = newPoint.tau + pdelta; // phase change
                    newPoint.q   = newPoint.q
                        + sddelta * rddelta * si * o.ccpx.real()
                            * oldPoint.p; // beam-width
                                          // change
                }
            }
        }
    } else {
        RunError(errState, BHC_ERR_BOUNDARY_CONDITION_TYPE);
        newPoint.Amp   = NAN;
        newPoint.Phase = NAN;
    }

    // printf("Reflection amp changed from to %g %g\n", oldPoint.Amp, newPoint.Amp);
}

} // namespace bhc
