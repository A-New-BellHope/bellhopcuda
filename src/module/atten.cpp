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
#include "atten.hpp"

namespace bhc { namespace module {

/**
 * Francois Garrison formulas for attenuation
 * Based on a Matlab version by D. Jackson APL-UW
 *
 * mbp Feb. 2019
 * Verified using F-G Table IV
 *
 * alpha = attenuation   (dB/km)
 * f     = frequency     (kHz)
 * T     = temperature   (deg C)
 * S     = salinity      (psu)
 * pH    = 7 for neutral water
 * z_bar = depth         (m)
 *
 *     Returns
 *        alpha = volume attenuation in dB/km
 */
inline real Franc_Garr(real f, const AttenInfo *atten)
{
    real c, a1, a2, a3, p1, p2, p3, f1, f2;

    c = FL(1412.0) + FL(3.21) * atten->t + FL(1.19) * atten->Salinity
        + FL(0.0167) * atten->z_bar;

    // Boric acid contribution
    a1 = FL(8.86) / c * STD::pow(FL(10.0), FL(0.78) * atten->pH - FL(5.0));
    p1 = FL(1.0);
    f1 = FL(2.8) * STD::sqrt(atten->Salinity / 35)
        * STD::pow(FL(1.0), FL(4.0) - FL(1245.0) / (atten->t + FL(273.0)));

    // Magnesium sulfate contribution
    a2 = FL(21.44) * atten->Salinity / c * (FL(1.0) + FL(0.025) * atten->t);
    p2 = FL(1.0) - RL(1.37e-4) * atten->z_bar + RL(6.2e-9) * SQ(atten->z_bar);
    f2 = FL(8.17) * STD::pow(FL(10.0), FL(8.0) - FL(1990.0) / (atten->t + 273))
        / (FL(1.0) + FL(0.0018) * (atten->Salinity - FL(35.0)));

    // Viscosity
    p3 = FL(1.0) - RL(3.83e-5) * atten->z_bar + RL(4.9e-10) * SQ(atten->z_bar);
    if(atten->t < FL(20.0)) {
        a3 = RL(4.937e-4) - RL(2.59e-5) * atten->t + RL(9.11e-7) * SQ(atten->t)
            - RL(1.5e-8) * CUBE(atten->t);
    } else {
        a3 = RL(3.964e-4) - RL(1.146e-5) * atten->t + RL(1.45e-7) * SQ(atten->t)
            - RL(6.5e-10) * CUBE(atten->t);
    }

    return a1 * p1 * (f1 * SQ(f)) / (SQ(f1) + SQ(f))
        + a2 * p2 * (f2 * SQ(f)) / (SQ(f2) + SQ(f)) + a3 * p3 * SQ(f);
}

/**
 * Converts real wave speed and attenuation to a single
 * complex wave speed (with positive imaginary part)
 *
 * AttenUnit
 * 6 Cases:    N Nepers/meter
 *             M dB/meter      (M for Meters)
 *             m dB/meter with a power law
 *             F dB/m-kHz      (F for frequency dependent)
 *             W dB/wavelength (W for Wavelength)
 *             Q Q
 *             L Loss parameter
 *
 * second letter adds volume attenuation according to standard laws:
 *             T for Thorp
 *             F for Francois Garrison
 *             B for biological
 *
 * freq is the current frequency
 * freq0 is the reference frequency for which the dB/meter was specified (used only for
 * 'm')
 *
 * Returns
 * c     real      part of sound speed
 * alpha imaginary part of sound speed
 */
template<bool O3D> cpx crci(
    const bhcParams<O3D> &params, real z, real c, real alpha, const char (&AttenUnit)[2])
{
    real f2, omega, alphaT, Thorp, a, fg;
    cpx ret;

    real freq, freq0;
    freq = freq0                = params.freqinfo->freq0;
    const AttenInfo *atten      = params.atten;
    PrintFileEmu &PRTFile       = GetInternal(params)->PRTFile;
    constexpr real betaPowerLaw = FL(1.0);

    omega = FL(2.0) * REAL_PI * freq;

    // Convert to Nepers/m
    alphaT = FL(0.0);
    switch(AttenUnit[0]) {
    case 'N': alphaT = alpha; break;
    case 'M': // dB/m
        alphaT = alpha / RL(8.6858896);
        break;
    case 'm': // dB/m with power law
        alphaT = alpha / RL(8.6858896);
        if(freq < params.fT) { // frequency raised to the power betaPowerLaw
            alphaT *= STD::pow(freq / freq0, betaPowerLaw);
        } else { // linear in frequency
            alphaT *= (freq / freq0)
                * STD::pow(params.fT / freq0, betaPowerLaw - FL(1.0));
        }
        break;
    case 'F': // dB/(m kHz)
        alphaT = alpha * freq / RL(8685.8896);
        break;
    case 'W': // dB/wavelength
        if(c != FL(0.0)) alphaT = alpha * freq / (RL(8.6858896) * c);
        // The following lines give f^1.25 frequency dependence
        // real FAC = STD::sqrt(STD::sqrt(freq / FL(50.0)));
        // if(c != FL(0.0)) alphaT = FAC * alpha * freq / (RL(8.6858896) * c);
        break;
    case 'Q': // Quality factor
        if(c * alpha != FL(0.0)) alphaT = omega / (FL(2.0) * c * alpha);
        break;
    case 'L': // loss parameter
        if(c != FL(0.0)) alphaT = alpha * omega / c;
        break;
    }

    // added volume attenuation

    switch(AttenUnit[1]) {
    case 'T':
        f2 = SQ(freq / FL(1000.0));

        // Original formula from Thorp 1967
        // Thorp = FL(40.0) * f2 / (FL(4100.0) + f2) + FL(0.1) * f2 / (FL(1.0) + f2);   //
        // dB/kyard Thorp /= RL(914.4);                 // dB / m Thorp /= RL(8.6858896);
        // // Nepers / m

        // Updated formula from JKPS Eq. 1.34
        Thorp = FL(3.3e-3) + FL(0.11) * f2 / (FL(1.0) + f2)
            + FL(44.0) * f2 / (FL(4100.0) + f2) + FL(3e-4) * f2; // dB/km
        Thorp /= FL(8685.8896);                                  // Nepers / m

        alphaT += Thorp;
        break;
    case 'F':                                      // Francois-Garrison
        fg = Franc_Garr(freq / FL(1000.0), atten); // dB/km
        fg /= FL(8685.8896);                       // Nepers / m
        alphaT += fg;
        break;
    case 'B': // biological attenuation per Orest Diachok
        for(int32_t iBio = 0; iBio < atten->NBioLayers; ++iBio) {
            if(z >= atten->bio[iBio].z1 && z <= atten->bio[iBio].z2) {
                a = atten->bio[iBio].a0
                    / (SQ(1.0 - SQ(atten->bio[iBio].f0) / SQ(freq))
                       + 1.0 / SQ(atten->bio[iBio].q)); // dB/km
                a /= FL(8685.8896);                     // Nepers / m
                alphaT += a;
            }
        }
        break;
    }

    // Convert Nepers/m to equivalent imaginary sound speed
    alphaT *= SQ(c) / omega;
    ret = cpx(c, alphaT);

    if(alphaT > c) {
        PRTFile << "Complex sound speed: " << ret << "\n";
        PRTFile << "Usually this means you have an attenuation that is way too high\n";
        EXTWARN("attenuation: crci: The complex sound speed has an imaginary part > "
                "real part");
    }

    return ret;
}

#if BHC_ENABLE_2D
template cpx crci<false>(
    const bhcParams<false> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
#endif
#if BHC_ENABLE_NX2D || BHC_ENABLE_3D
template cpx crci<true>(
    const bhcParams<true> &params, real z, real c, real alpha,
    const char (&AttenUnit)[2]);
#endif

}} // namespace bhc::module
