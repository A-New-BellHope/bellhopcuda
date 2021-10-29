#include "common.hpp"

struct bioStructure {
    real z1, z2, f0, q, a0;
};

constexpr int32_t MaxBioLayers = 200;

struct AttenInfo {
    int32_t NBioLayers;
    bioStructure bio[MaxBioLayers];
    real t, Salinity, pH, z_bar, fg;
};

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
real Franc_Garr(real f, const AttenInfo *atten){
    real c, a1, a2, a3, p1, p2, p3, f1, f2;
    
    c = RC(1412.0) + RC(3.21) * atten->t + RC(1.19) * atten->Salinity + RC(0.0167) * atten->z_bar;

    // Boric acid contribution
    a1 = RC(8.86) / c * STD::pow(RC(10.0), RC(0.78) * pH - RC(5.0));
    p1 = RC(1.0);
    f1 = RC(2.8) * STD::sqrt(atten->Salinity / 35) * STD::pow(RC(1.0), RC(4.0) - RC(1245.0) / (t + RC(273.0)));

    // Magnesium sulfate contribution
    a2 = RC(21.44) * atten->Salinity / c * (RC(1.0) + RC(0.025) * atten->t);
    p2 = RC(1.0) - RC(1.37e-4) * atten->z_bar + RC(6.2e-9) * SQ(atten->z_bar);
    f2 = RC(8.17) * STD::pow(RC(10.0), RC(8.0) - RC(1990.0) / (atten->t + 273)) / (RC(1.0) + RC(0.0018) * (atten->Salinity - RC(35.0)));

    // Viscosity
    p3 = RC(1.0) - RC(3.83e-5) * atten->z_bar + RC(4.9e-10) * SQ(atten->z_bar);
    if(t < RC(20.0)){
        a3 = RC(4.937e-4) - RC(2.59e-5) * atten->t + RC(9.11e-7) * SQ(atten->t) - RC(1.5e-8) * CUBE(atten->t);
    }else{
        a3 = RC(3.964e-4) RC(-1.146e-5) * atten->t + RC(1.45e-7) * SQ(atten->t) - RC(6.5e-10) * CUBE(atten->t);
    }
    
    return a1 * p1 * (f1 * SQ(f)) / (SQ(f1) + SQ(f)) 
         + a2 * p2 * (f2 * SQ(f)) / (SQ(f2) + SQ(f))
         + a3 * p3 * SQ(f)
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
 * freq0 is the reference frequency for which the dB/meter was specified (used only for 'm')
 * 
 * Returns
 * c     real      part of sound speed
 * alpha imaginary part of sound speed
 */
cpx crci(real z, real c, real alpha, real freq, real freq0, char (&AttenUnit)[2],
    real beta, real fT, const AttenInfo *atten, std::ofstream &PRTFile)
{
    real f2, omega, alphaT, Thorp, a, fg;
    cpx ret;
    
    omega = RC(2.0) * M_PI * freq;
    
    // Convert to Nepers/m
    alphaT = RC(0.0);
    switch(AttenUnit[0]){
    case 'N':
        alphaT = alpha; break;
    case 'M': // dB/m
        alphaT = alpha / RC(8.6858896); break;
    case 'm': // dB/m with power law
        alphaT = alpha / RC(8.6858896);
        if(freq < fT){ // frequency raised to the power beta
            alphaT *= STD::pow(freq / freq0, beta);
        }else{ // linear in frequency
            alphaT *= (freq / freq0) * STD::pow(fT / freq0, beta - RC(1.0));
        }
        break;
    case 'F': // dB/(m kHz)
        alphaT = alpha * freq / RC(8685.8896); break;
    case 'W': // dB/wavelength
        if(c != RC(0.0)) alphaT = alpha * freq / (RC(8.6858896) * c);
        // The following lines give f^1.25 frequency dependence
        // real FAC = STD::sqrt(STD::sqrt(freq / RC(50.0)));
        // if(c != RC(0.0)) alphaT = FAC * alpha * freq / (RC(8.6858896) * c);
        break;
    case 'Q': // Quality factor
        if(c * alpha != RC(0.0)) alphaT = omega / (RC(2.0) * c * alpha);
        break;
    case 'L': // loss parameter
        if(c != RC(0.0)) alphaT = alpha * omega / c;
        break;
    }
    
    // added volume attenuation
    
    switch(AttenUnit[1]){
    case 'T':
        f2 = SQ(freq / RC(1000.0));
        
        // Original formula from Thorp 1967
        // Thorp = RC(40.0) * f2 / (RC(4100.0) + f2) + RC(0.1) * f2 / (RC(1.0) + f2);   // dB/kyard
        // Thorp /= RC(914.4);                 // dB / m
        // Thorp /= RC(8.6858896);             // Nepers / m
        
        // Updated formula from JKPS Eq. 1.34
        Thorp = RC(3.3e-3) + RC(0.11) * f2 / (RC(1.0) + f2) + 
            RC(44.0) * f2 / (RC(4100.0) + f2) + RC(3e-4) * f2; // dB/km
        Thorp /= RC(8685.8896); // Nepers / m
 
        alphaT += Thorp;
        break;
    case 'F': // Francois-Garrison
        fg = Franc_Garr(freq / RC(1000.0), atten); // dB/km
        fg /= RC(8685.8896); // Nepers / m
        alphaT += fg;
        break;
    case 'B': // biological attenuation per Orest Diachok
        for(int32_t iBio=0; iBio<atten->NBioLayers; ++iBio){
            if(z >= atten->bio[iBio].z1 && z <= atten->bio[iBio].z2){
                a = atten->bio[iBio].a0 / (SQ(1.0 - SQ(atten->bio[iBio].f0) / SQ(freq))
                    + 1.0 / SQ(atten->bio[iBio].q)); // dB/km
                a /= RC(8685.8896); // Nepers / m
                alphaT += a;
            }
        }
        break;
    }
    
    // Convert Nepers/m to equivalent imaginary sound speed 
    alphaT *= SQ(c) / omega;
    ret = cpx(c, alphaT);
    
    if(alphaT > c){
        PRTFile << "Complex sound speed: " << ret << "\n";
        PRTFile << "Usually this means you have an attenuation that is way too high\n";
        std::cout << "attenuation: crci: The complex sound speed has an imaginary part > real part\n";
    }
    
    return ret;
}
