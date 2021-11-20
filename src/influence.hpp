#pragma once
#include "common.hpp"
#include "atomics.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"
#include "step.hpp"

struct InfluenceRayInfo {
    // LP: Variables carried over between iterations.
    real phase, qOld, q0;
    vec2 x;
    bool lastValid;
    real zn0, rn0;
    int32_t kmah;
    int32_t ir;
    cpx gamma;
    // LP: Constants.
    cpx epsilon; // beam constant
    real freq0, omega;
    real RadMax;
    int32_t iBeamWindow2, NRz_per_range;
    real SrcDeclAngle; // take-off angle in degrees
    real Dalpha; // angular spacing
    real Ratio1; // scale factor (point source vs. line source)
};

/**
 * Calculates a smoothing function based on the h0 hermite cubic
 * x is the point where the function is to be evaluated
 * returns:
 * [ 0, x1 ] = 1
 * [x1, x2 ] = cubic taper from 1 to 0
 * [x2, inf] = 0
 */
HOST_DEVICE inline real Hermite(real x, real x1, real x2)
{
    real Ax = STD::abs(x);
    if(Ax <= x1){
        return RC(1.0);
    }else if(Ax >= x2){
        return RC(0.0);
    }else{
        real u = (Ax - x1) / (x2 - x1);
        return (RC(1.0) + RC(2.0) * u) * SQ(RC(1.0) - u);
    }
    // ret /= (RC(0.5) * (x1 + x2));
}

/**
 * Scale the pressure field
 * 
 * r: ranges (LP: [Nr])
 * Dalpha: angular spacing between rays
 * freq: source frequency
 * c: nominal sound speed (LP: [NRz][Nr])
 * u: Pressure field
 */
HOST_DEVICE inline void ScalePressure(real Dalpha, real c, real *r, 
    cpx *u, int32_t NRz, int32_t Nr, char (&RunType)[5], real freq)
{
    // Compute scale factor for field
    real cnst;
    if(RunType[1] == 'C' || RunType[1] == 'R'){
        // Cerveny Gaussian beams in Cartesian or Ray-centered coordinates
        cnst = -Dalpha * STD::sqrt(freq) / c;
    }else{
        cnst = RC(-1.0);
    }
    
    // For incoherent run, convert intensity to pressure
    if(RunType[0] != 'C'){
        for(int32_t irz=0; irz<NRz; ++irz){
            for(int32_t ir=0; ir<Nr; ++ir){
                u[irz*Nr+ir] = cpx(STD::sqrt(u[irz*Nr+ir].real()), RC(0.0));
            }
        }
    }
    
    // scale and/or incorporate cylindrical spreading
    for(int32_t ir=0; ir<Nr; ++ir){
        real factor;
        if(RunType[3] == 'X'){ // line source
            factor = RC(-4.0) * STD::sqrt(M_PI) * cnst;
        }else{ // point source
            if(r[ir] == RC(0.0)){
                factor = RC(0.0); // avoid /0 at origin, return pressure = 0
            }else{
                factor = cnst / STD::sqrt(STD::abs(r[ir]));
            }
        }
        for(int32_t irz=0; irz<NRz; ++irz){
            u[irz*Nr+ir] *= factor;
        }
    }
}

/**
 * Checks for a branch cut crossing and updates kmah accordingly
 */
HOST_DEVICE inline void BranchCut(const cpx &q1C, const cpx &q2C,
    char (&BeamType)[4], int32_t &kmah)
{
    real q1, q2;
    if(BeamType[1] == 'W'){ // WKBeams
        q1 = q1C.real();
        q2 = q2C.real();
    }else{
        if(q2C.real() >= RC(0.0)) return;
        q1 = q1C.imag();
        q2 = q2C.imag();
    }
    if( (q1 < RC(0.0) && q2 >= RC(0.0)) || 
        (q1 > RC(0.0) && q2 <= RC(0.0)) ) kmah = -kmah;
}

HOST_DEVICE inline void ApplyContribution(real cnst, real w,
    real omega, cpx delay, real phaseInt, real SrcDeclAngle, real RcvrDeclAngle,
    cpx *u, const BeamStructure *Beam)
{
    switch(Beam->RunType[0]){
    case 'E':
        // eigenrays
        print("Eigenrays not yet supported\n");
        bail();
        //WriteRay2D(SrcDeclAngle, is);
        break;
    case 'A':
    case 'a':
        // arrivals
        print("Arrivals not yet supported\n");
        bail();
        //AddArr(omega, iz, ir, Amp, phaseInt, delay, SrcDeclAngle, RcvrDeclAngle,
        //    point1.NumTopBnc, point1.NumBotBnc);
        break;
    case 'C':
        // coherent TL
        AtomicAddCpx(u, cnst * w * STD::exp(-J * (omega * delay - phaseInt)));
        // omega * SQ(n) / (RC(2.0) * SQ(point1.c) * delay)))) // curvature correction
        break;
    default:
        // incoherent/semicoherent TL
        real v = STD::pow(cnst * STD::exp((omega * delay).imag()), RC(2.0) * w);
        if(Beam->Type[0] == 'B'){
            // Gaussian beam
            v *= STD::sqrt(RC(2.0) * M_PI);
        }
        AtomicAddCpx(u, cpx(v, RC(0.0)));
    }
}

/**
 * Picks the optimum value for epsilon
 * 
 * omega: angular frequency
 * c: sound speed
 * gradc: gradient
 * alpha: angular spacing for ray fan
 * rLoop: loop range
 * EpsMultiplier: multiplier [LP: :( ]
 */
HOST_DEVICE inline cpx PickEpsilon(char BeamType0, char BeamType1, real omega, 
    real c, vec2 gradc, real alpha, real Dalpha, real rLoop, real EpsMultiplier)
{
    real halfwidth, cz;
    cpx epsilonOpt;
    //const char *tag;
    switch(BeamType0){
    case 'C':
    case 'R':
        //tag = "Paraxial beams";
        switch(BeamType1){
        case 'F':
            //tag = "Space filling beams";
            halfwidth = RC(2.0) / ((omega / c) * Dalpha);
            epsilonOpt = J * RC(0.5) * omega * SQ(halfwidth);
            break;
        case 'M':
            //tag = "Minimum width beams";
            halfwidth = STD::sqrt(RC(2.0) * c * RC(1000.0) * rLoop / omega);
            epsilonOpt = J * RC(0.5) * omega * SQ(halfwidth);
            break;
        case 'W':
            //tag = "WKB beams";
            halfwidth = REAL_MAX;
            cz = gradc.y;
            if(cz == RC(0.0)){
                epsilonOpt = RC(1e10);
            }else{
                epsilonOpt = (-STD::sin(alpha) / STD::cos(SQ(alpha))) * c * c / cz;
            }
            break;
        default:
            print("Invalid BeamType[1]: %c\n", BeamType1);
            bail();
        }
        break;
    case 'G':
    case 'g':
        //tag = "Geometric hat beams";
        halfwidth = RC(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * RC(0.5) * omega * SQ(halfwidth);
        break;
    case 'B':
        //tag = "Geometric Gaussian beams";
        halfwidth = RC(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * RC(0.5) * omega * SQ(halfwidth);
        break;
    case 'b':
        printf("bellhopcuda: Geo Gaussian beams in ray-cent. coords. not "
            "implemented in BELLHOP (and therefore not in bellhopcuda)\n");
        bail();
        break;
    case 'S':
        //tag = "Simple Gaussian beams";
        halfwidth = RC(2.0) / ((omega / c) * Dalpha);
        epsilonOpt = J * RC(0.5) * omega * SQ(halfwidth);
        break;
    default:
        print("Invalid BeamType[0]: %c\n", BeamType0);
        bail();
    }
    
    /*
    // On first call write info to prt file
    static bool INIFlag = true;
    if(INIFlag){
        PRTFile << "\n" << tag << "\n";
        PRTFile << "halfwidth  = " << halfwidth << "\n";
        PRTFile << "epsilonOpt = " << epsilonOpt << "\n";
        PRTFile << "EpsMult    = " << EpsMultiplier << "\n\n";
        INIFlag = false;
    }
    */
    
    return EpsMultiplier * epsilonOpt;
}

// LP: Helper functions.

/**
 * phase shifts at caustics
 */
HOST_DEVICE inline void IncPhaseIfCaustic(InfluenceRayInfo &inflray, real q)
{
    if(q < RC(0.0) && inflray.qOld >= RC(0.0) || q > RC(0.0) && inflray.qOld <= RC(0.0))
        inflray.phase += M_PI * RC(0.5);
}

/**
 * phase shifts at caustics
 * LP: BUG: point.phase is discarded if the condition is met.
 */
HOST_DEVICE inline real BuggyFinalPhase(const ray2DPt &point, 
    const InfluenceRayInfo &inflray, real q)
{
    real phaseInt = point.Phase + inflray.phase;
    if(q <= RC(0.0) && inflray.qOld > RC(0.0) || q >= RC(0.0) && inflray.qOld < RC(0.0))
        phaseInt = inflray.phase + M_PI * RC(0.5);
    return phaseInt;
}

HOST_DEVICE inline int32_t RToIR(real r, const Position *Pos)
{
    // mbp: should be ", -1);"? [LP: for computation of ir1 / irA]
    return STD::max(STD::min((int)((r - Pos->Rr[0]) / Pos->Delta_r), Pos->NRr-1), 0);
}

HOST_DEVICE inline void Compute_N_R_IR(real &n, real &r, int32_t &ir,
    real xx, real xy, real zn, real rn, real zR, const Position *Pos)
{
    n = (zR - xy) / zn;
    r = xx + n * rn;
    // following assumes uniform spacing in Pos->r
    ir = RToIR(r, Pos); // index of receiver
}

/**
 * detect and skip duplicate points (happens at boundary reflection)
 */
HOST_DEVICE inline bool IsDuplicatePoint(const ray2DPt &point0, const ray2DPt &point1)
{
    return STD::abs(point1.x.x - point0.x.x) < RC(1.0e3) * spacing(point1.x.x);
}

HOST_DEVICE inline void Compute_eps_pB_qB(cpx &eps, cpx &pB, cpx &qB,
    const ray2DPt &point, const InfluenceRayInfo &inflray, const BeamStructure *Beam)
{
    if(Beam->Type[1] == 'C'){
        eps = J * STD::abs(point.q.x / point.q.y);
    }else{
        eps = inflray.epsilon;
    }
    pB = point.p.x + eps * point.p.y;
    qB = point.q.x + eps * point.q.y;
}

// 

HOST_DEVICE inline void Init_Influence(InfluenceRayInfo &inflray,
    const ray2DPt &point0, real alpha, vec2 gradc,
    const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam)
{
    if(Beam->Type == 'R' || Beam->Type == 'C'){
        inflray.epsilon = PickEpsilon(Beam->Type[0], Beam->Type[1], inflray.omega,
            point0.c, gradc, alpha, Angles->Dalpha, Beam->rLoop, Beam->epsMultiplier);
    }
    inflray.freq0 = freqinfo->freq0;
    inflray.omega = RC(2.0) * M_PI * inflray.freq0;
    // LP: The 5x version is changed to 50x on both codepaths before it is used.
    // inflray.RadMax = RC(5.0) * ccpx.real() / freqinfo->freq0; // 5 wavelength max radius
    inflray.RadMax = RC(50.0) * point0.c / freqinfo->freq0; // 50 wavelength max radius
    inflray.iBeamWindow2 = SQ(Beam->iBeamWindow);
    inflray.NRz_per_range = Compute_NRz_per_range(Pos, Beam);
    inflray.SrcDeclAngle = RadDeg * alpha;
    inflray.Dalpha = Angles->Dalpha;
    if(Beam->RunType[3] == 'R'){
        inflray.Ratio1 = STD::sqrt(STD::abs(STD::cos(alpha))); // point source
        // LP: In InfluenceSGB, this was sqrt(cos(alpha)) without the abs,
        // but Ratio1 is real, so there will be an error for backward rays.
        // I don't think there's any real reason this influence function should
        // not be capable of handling backward rays, so this check could just
        // be removed, but we must match all the quirks of BELLHOP.
        if(STD::cos(alpha) < RC(0.0) && Beam->Type[0] == 'S'){
            printf("Error: Original SGB implementation has spurious error for backward rays\n");
            bail();
        }
    }else{
        inflray.Ratio1 = RC(1.0); // line source
    }
    
    // LP: For GeoHat and Gaussian
    inflray.q0 = point0.c / inflray.Dalpha; // Reference for J = q0 / q
    
    // LP: For all except Cerveny
    inflray.phase = RC(0.0);
    
    if(Beam->Type[0] == 'S'){
        // LP: For SGB
        inflray.qOld = RC(1.0);
    }else{
        // LP: For GeoHat and Gaussian
        inflray.qOld = point0.q.x; // used to track KMAH index
    }
    
    // LP: For Cerveny
    inflray.kmah = 1;
    
    // LP: For RayCen types
    if(Beam->Type[0] == 'g'){
        // LP: For GeoHat
        inflray.zn = -point0.t.x * point0.c;
        inflray.rn =  point0.t.y * point0.c;
        inflray.x = point0.x;
        inflray.lastValid = STD::abs(inflray.zn) >= 1e-6;
    }else{
        // LP: For Cerveny
        // mbp: This logic means that the first step along the ray is skipped
        // which is a problem if deltas is very large, e.g. isospeed problems
        // I fixed this in InfluenceGeoHatRayCen
        inflray.zn = RC(1.0);
        inflray.rn = RC(0.0);
        inflray.x = vec2(RC(0.0), RC(0.0));
        inflray.lastValid = false;
    }
}

// LP: Main functions.

/**
 * Paraxial (Cerveny-style) beams in ray-centered coordinates
 */
HOST_DEVICE inline void Step_InfluenceCervenyRayCen(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpx *u,
    const BdryType *Bdry, const Position *Pos, const BeamStructure *Beam)
{
    cpx eps0, eps1, qB0, qB1, gamma0, gamma1;
    real zn, rn, zR;
    // need to add logic related to NRz_per_range
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB(eps1, pB1, qB1, point1, inflray, Beam);
    gamma0 = pB0 / qB0;
    gamma1 = pB1 / qB1;
    
    // ray normal based on tangent with c(s) scaling
    zn = -point1.t.x * point1.c;
    rn =  point1.t.y * point1.c;
    // If normal parallel to TL-line, skip to next step on ray
    // LP: Possible BUG: This is the same as abs(zn) <= 0 -- there are no other
    // numbers less than this. Definitely not equivalent to something like
    // abs(zn) < RC(1e-7). The function was originally FORTRAN's TINY.
    if(STD::abs(zn) < REAL_MINPOS) return;
    
    // detect and skip duplicate points (happens at boundary reflection)
    if(IsDuplicatePoint(point0, point1)){
        inflray.lastValid = true;
        inflray.x = point1.x;
        inflray.zn = zn;
        inflray.rn = rn;
        return;
    }
    
    // compute KMAH index
    // Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
    int32_t old_kmah = inflray.kmah;
    BranchCut(qB0, qB1, Beam->Type, inflray.kmah);
    
    for(int32_t iz=0; iz<inflray.NRz_per_range; ++iz){
        zR = Pos->Rz[iz];
        
        for(int32_t image=1; image <= Beam->Nimage; ++image){
            real tempz = point1.x.y;
            real temprn = rn;
            // LP: BUG: The original code inverted rn for beam 2 and inverted it
            // again for beam 3, meaning that if Beam->Nimage is set to 2, rn will
            // be inverted for every other step.
            if(image == 1){ // True beam
                (void)0;
            }else if(image == 2){ // Surface-reflected beam
                tempz = RC(2.0) * Bdry->Top.hs.Depth - tempz;
                temprn = -temprn;
            }else if(image == 3){ // Bottom-reflected beam
                tempz = RC(2.0) * Bdry->Bot.hs.Depth - tempz;
            }else{
                printf("Invalid Beam->Nimage %d\n", Beam->Nimage);
                bail();
            }
            real nA, rA, nB, rB;
            int32_t ir1, ir2; // LP: mbp switches from A/B naming to 1/2 here.
            Compute_N_R_IR(nB, rB, ir2, point1.x.x, tempz, zn, temprn, zR, Pos);
            Compute_N_R_IR(nA, rA, ir1, inflray.x.x, inflray.x.y, inflray.zn, 
                inflray.rn, zR, Pos);
            
            if(inflray.lastValid && ir1 < ir2){
                for(int32_t ir=ir1+1; ir<=ir2; ++ir){
                    real w, n, nSq, c;
                    cpx q, gamma, tau, contri;
                    w     = (Pos->Rr[ir] - rA) / (rB - rA);
                    q     = qB0    + w * (qB1 - qB0);
                    gamma = gamma0 + w * (gamma1 - gamma0);
                    n     = nA     + w * (nB - nA);
                    nSq   = SQ(n);
                    if(gamma.imag() > 0){
                        print("Unbounded beam\n");
                        continue;
                    }
                    
                    if(RC(-0.5) * inflray.omega * gamma.imag() * nSq < inflray.iBeamWindow2){ // Within beam window?
                        c   = point0.c;
                        tau = point0.tau + w * (point1.tau - point0.tau);
                        contri = inflray.Ratio1 * point1.Amp * STD::sqrt(c * STD::abs(eps1) / q) *
                            STD::exp(-J * (inflray.omega * (tau + RC(0.5) * gamma * nSq) - point1.phase));
                        
                        cpx P_n = -J * inflray.omega * gamma * n * contri;
                        cpx P_s = -J * inflray.omega / c         * contri;
                        switch(Beam->Component){
                        case 'P': // pressure
                        case 'V': // vertical component
                            contri = c * glm::dot(vec2(P_n, P_s), point1.t);
                            break;
                        case 'H': // horizontal component
                            contri = c * (-P_n * point1.t.y + P_s * point1.t.x);
                            break;
                        }
                        
                        int32_t kmah = old_kmah;
                        BranchCut(qB0, q, Beam->Type, kmah); // Get correct branch of STD::sqrt
                        
                        if(kmah < 0)   contri = -contri;
                        if(image == 2) contri = -contri;
                        
                        if(Beam->RunType[0] == 'I' || Beam->RunType[0] == 'S'){ // Incoherent or Semi-coherent TL
                            contri = contri * contri.conj();
                        }
                        
                        AtomicAddCpx(&u[iz*Pos->NRr+ir], 
                            Hermite(n, inflray.RadMax, RC(2.0) * inflray.RadMax) * contri);
                    }
                }
            }
        }
    }
    
    inflray.lastValid = true;
    inflray.x = point1.x;
    inflray.zn = zn;
    inflray.rn = rn;
}

/**
 * Paraxial (Cerveny-style) beams in ray-centered coordinates
 */
HOST_DEVICE inline void Step_InfluenceCervenyCart(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpx *u, 
    const BdryType *Bdry, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr, 
    const Position *Pos, const BeamStructure *Beam)
{
    cpx eps0, eps1, pB0, pB1, qB0, qB1, gamma0, gamma1;
    real zn, rn, zR;
    // need to add logic related to NRz_per_range
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    Compute_eps_pB_qB(eps0, pB0, qB0, point0, inflray, Beam);
    Compute_eps_pB_qB(eps1, pB1, qB1, point1, inflray, Beam);
    
    // Form gamma and KMAH index
    // Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType
    vec2 rayt = point1.c * point1.t; // unit tangent
    vec2 rayn = vec2(rayt.y, -rayt.x); // unit normal
    
    cpx ccpx; vec2 gradc; real crr, crz, czz, rho;
    EvaluateSSP(point1.x, ccpx, gradc, crr, crz, czz, rho, inflray.freq0, ssp, iSegz, iSegr);
    
    real csq = SQ(ccpx.real());
    real cS = glm::dot(gradc, rayt);
    real cN = glm::dot(gradc, rayn);
    
    real Tr = rayt.x;
    real Tz = rayt.y;
    
    gamma1 = RC(0.0);
    if(qB1 != RC(0.0)) gamma = RC(0.5) * (pB1 / qB1 * SQ(Tr) + 
        RC(2.0) * cN / csq * Tz * Tr - cS / csq * SQ(Tz));
    gamma0 = inflray.gamma;
    inflray.gamma = gamma1;
        
    int32_t old_kmah = inflray.kmah;
    BranchCut(qB0, qB1, Beam->Type, inflray.kmah);
    
    if(is < 2) return; // LP: Skips the first valid pair, which here would be is=1.
    if(point1.x.x > Pos->Rr[Pos->NRr-1]) return;
    real rA = point0.x.x;
    real rB = point1.x.x;
    if(IsDuplicatePoint(point0, point1)) return; // don't process duplicate points
    
    // Compute upper index on rcvr line
    // Assumes r is a vector of equally spaced points
    int32_t irA = RToIR(rA, Pos);
    int32_t irB = RToIR(rB, Pos);
    
    if(irA >= irB) return;
    
    for(int32_t ir=ir1+1; ir<=ir2; ++ir){
        real w, c;
        vec2 x, rayt;
        cpx q, tau, gamma, cnst;                
        w     = (Pos->Rr[ir] - rA) / (rB - rA);
        x     = point0.x   + w * (point1.x   - point0.x);
        rayt  = point0.t   + w * (point1.t   - point0.t);
        c     = point0.c   + w * (point1.c   - point0.c);
        q     = qB0        + w * (qB1        - qB0);
        tau   = point0.tau + w * (point1.tau - point0.tau);
        gamma = gamma0     + w * (gamma1     - gamma0);
        
        if(gamma.imag() > 0){
            print("Unbounded beam\n");
            continue;
        }
        
        cnst = inflray.Ratio1 * STD::sqrt(c * STD::abs(eps0) / q);
        
        // Get correct branch of STD::sqrt
        int32_t kmah = old_kmah;
        BranchCut(qB0, q, Beam->Type, kmah);
        if(kmah < 0) cnst = -cnst;
        
        for(int32_t iz=0; iz<inflray.NRz_per_range; ++iz){
            zR = Pos->Rz[iz];
            
            cpx contri = RC(0.0);
            for(int32_t image=1; image <= Beam->Nimage; ++image){
                real deltaz, Polarity;
                if(image == 1){ // True beam
                    deltaz = zR - x.y;
                    Polarity = RC(1.0);
                }else if(image == 2){ // Surface reflected beam
                    deltaz = -zR + RC(2.0) * Bdry->Top.hs.Depth - x.y;
                    Polarity = RC(-1.0);
                }else if(image == 3){ // Bottom  reflected beam
                    deltaz = -zR + RC(2.0) * Bdry->Bot.hs.Depth - x.y;
                    Polarity = RC(1.0); // assumes rigid bottom
                }else{
                    printf("Invalid Beam->Nimage %d\n", Beam->Nimage);
                    bail();
                }
                if(inflray.omega * gamma.imag() * SQ(deltaz) < inflray.iBeamWindow2)
                    contri += Polarity * point1.Amp * 
                        Hermite(deltaz, inflray.RadiusMax, RC(2.0) * inflray.RadiusMax) *
                        STD::exp(-J * (inflray.omega * (
                            tau + rayt.y * deltaz + gamma * SQ(deltaz)) - point1.Phase));
            }
            
            if(Beam->RunType[0] == 'C'){ // coherent
                contri = cnst * contri;
            }else if(Beam->RunType[0] == 'I' || Beam->RunType[0] == 'S'){
                contri = cnst * contri;
                contri = contri * contri.conj();
            }
            AtomicAddCpx(&u[iz*Pos->NRr+ir], contri);
        }
    }
}

/**
 * Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates
 */
HOST_DEVICE inline void Step_InfluenceGeoHatRayCen(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray,
    int32_t is, cpx *u, 
    const Position *Pos, const BeamStructure *Beam)
{
    real RcvrDeclAngle;
    real dq, q, w, n, l, delay, cnst, phaseInt;
    real zn, rn, nA, nB, rA, rB;
    int32_t irA, irB;
    cpx dtau;
    
    dq = point1.q.x - point0.q.x;
    dtau = point1.tau - point0.tau;
    
    // ray normal based on tangent with c(s) scaling
    zn = -point1.t.x * point1.c;
    rn =  point1.t.y * point1.c;
    if(STD::abs(zn) < RC(1e-10)) return;
    
    if(IsDuplicatePoint(point0, point1)){
        inflray.lastValid = true;
        inflray.x = point1.x;
        inflray.zn = zn;
        inflray.rn = rn;
        return;
    }
    
    RcvrDeclAngle = RadDeg * STD::atan2(point1.t.y, point1.t.x);
    
    // During reflection imag(q) is constant and adjacent normals cannot bracket
    // a segment of the TL line, so no special treatment is necessary
    
    real scaledAmp = inflray.Ratio1 * STD::sqrt(point1.c) * point1.Amp;
    
    for(int32_t iz=0; iz<inflray.NRz_per_range; ++iz){
        real zR = Pos->Rz[iz];
        
        if(!inflray.lastValid){
            nA = RC(1e10);
            rA = RC(1e10);
            irA = 0;
        }else{
            Compute_N_R_IR(nA, rA, irA, inflray.x.x, inflray.x.y, 
                inflray.zn, inflray.rn, zR, Pos);
        }
        Compute_N_R_IR(nB, rB, irB, point1.x.x, point1.x.y, zn, rn, zR, Pos);
        
        if(irA == irB) continue;
        
        q = point0.q.x;
        IncPhaseIfCaustic(inflray, q);
        inflray.qOld = q;
        
        // *** Compute contributions to bracketted receivers ***
        
        for(int32_t ir = STD::min(irA, irB) + 1; ir <= STD::max(irA, irB); ++ir){
            w = (Pos->Rr[ir] - rA) / (rB - rA);
            n = STD::abs(nA + w * (nB - nA));
            q = point0.q.x + w * dq; // interpolated amplitude
            l = STD::abs(q) / inflray.q0; // beam radius
            
            if(n < l){ // in beamwindow?
                delay = point0.tau + w * dtau;
                cnst  = scaledAmp / STD::sqrt(STD::abs(q));
                w     = (l - n) / l; // hat function: 1 on center, 0 on edge
                phaseInt = BuggyFinalPhase(point0, inflray, q);
                
                ApplyContribution(cnst, w, inflray.omega, 
                    delay, phaseInt, inflray.SrcDeclAngle, RcvrDeclAngle,
                    &u[iz*Pos->NRr+ir], Beam);
            }
        }
    }
    
    inflray.lastValid = true;
    inflray.x = point1.x;
    inflray.zn = zn;
    inflray.rn = rn;
}

/**
 * Geometric, hat-shaped or Gaussian beams in Cartesian coordintes
 *
 * u: complex pressure field
 */
HOST_DEVICE inline void Step_InfluenceGeoHatOrGaussianCart(
    bool isGaussian, 
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpx *u, 
    const Position *Pos, const BeamStructure *Beam)
{
    const int32_t BeamWindow = 4; // beam window: kills beams outside e**(-0.5 * ibwin**2 )
    
    vec2 x_ray, rayt, rayn, x_rcvr;
    real q0, RcvrDeclAngle;
    real rlen, RadiusMax, zMin, zMax, sigma, lambda, a, dqds;
    cpx dtauds;
    
    rA = point0.x.x;
    rB = point1.x.x;
    
    // what if never satistified?
    // what if there is a single receiver (ir = -1 possible)
    // LP: This implementation has been adjusted to handle these cases.
    int32_t ir = BinarySearchGEQ(Pos->Rr, Pos->NRr, 1, 0, STD::min(rA, rB));
    if(Pos->Rr[ir] < STD::min(rA, rB)) return; // not "bracketted"
    
    x_ray = point0.x;
    
    // compute normalized tangent (compute it because we need to measure the step length)
    rayt = point1.x - point0.x;
    rlen = glm::norm(tayt);
    if(rlen < RC(1.0e3) * spacing(point1.x.x)) return; // if duplicate point in ray, skip to next step along the ray
    rayt /= rlen;
    rayn = vec2(-rayt.y, rayt.x); // unit normal to ray
    RcvrDeclAngle = RadDeg * STD::atan2(rayt.y, rayt.x);
    
    dqds   = point1.q.x - point0.q.x;
    dtauds = point1.tau - point0.tau;
    
    q = point0.q.x;
    IncPhaseIfCaustic(inflray, q);
    inflray.qOld = q;
    
    sigma = STD::max(STD::abs(point0.q.x), STD::abs(point1.q.x) / 
        (inflray.q0 * STD::abs(rayt.x)); // beam radius projected onto vertical line
    if(isGaussian){
        // calculate beam width
        lambda    = point0.c / inflray.freq0;
        sigma     = STD::max(sigma, STD::min(RC(0.2) * inflray.freq0 * point1.tau, M_PI * lambda));
        RadiusMax = BeamWindow * sigma;
    }else{
        RadiusMax = sigma;
    }
    
    // depth limits of beam
    if(STD::abs(rayt.x) > RC(0.5)){ // shallow angle ray
        zmin = min(point0.x.y, point1.x.y) - RadiusMax;
        zmax = max(point0.x.y, point1.x.y) + RadiusMax;
    }else{
        zmin = -REAL_MAX;
        zmax =  REAL_MAX;
    }
    
    // compute beam influence for this segment of the ray
    for(; ir<Pos->NRr; ++ir){
        // is Rr[ir] contained in [rA, rB)? Then compute beam influence
        // LP: Because of the new setup and always incrementing regardless of
        // which direction the ray goes, we only have to check this side.
        if(Pos->Rr[ir] >= STD::max(rA, rB)) return;
        
        for(int32_t iz=0; iz<inflray.NRz_per_range; ++iz){
            if(Beam->RunType[4] == 'I'){
                x_rcvr = vec2(Pos->Rr[ir], Pos->Rz[ir]); // irregular grid
            }else{
                x_rcvr = vec2(Pos->Rr[ir], Pos->Rz[iz]); // rectilinear grid
            }
            if(x_rcvr.y < zmin || x_rcvr.x > zmax) continue;
            
            s     = glm::dot(x_rcvr - x_ray, rayt) / rlen; // proportional distance along ray
            n     = STD::abs(glm::dot(x_rcvr - x_ray, rayn)); // normal distance to ray
            q     = point0.q.x + s * dqds; // interpolated amplitude
            sigma = STD::abs(q / inflray.q0);
            real beamWCompare;
            if(isGaussian){
                sigma = STD::max(sigma, STD::min(RC(0.2) * inflray.freq0 * point1.tau, M_PI * lambda)); // min pi * lambda, unless near
                beamWCompare = BeamWindow * sigma;
            }else{
                RadiusMax = sigma;
                beamWCompare = RadiusMax;
            }
            
            if(n < beamWCompare){ // Within beam window?
                delay = point0.tau + s * dtauds; // interpolated delay
                cnst  = inflray.Ratio1 * STD::sqrt(point1.c / STD::abs(q)) * point1.Amp;
                if(isGaussian){
                    // sqrt( 2 * pi ) represents a sum of Gaussians in free space
                    cnst /= STD::sqrt(RC(2.0) * M_PI);
                    a = STD::abs(inflray.q0 / q);
                    w = STD::exp(RC(-0.5) * SQ(n / sigma)) / (sigma * a); // Gaussian decay
                }else{
                    w = (RadiusMax - n) / RadiusMax; // hat function: 1 on center, 0 on edge
                }
                phaseInt = BuggyFinalPhase((isGaussian ? point1 : point0), inflray, q);
                
                ApplyContribution(cnst, w, inflray.omega, 
                    delay, phaseInt, inflray.SrcDeclAngle, RcvrDeclAngle,
                    &u[iz*Pos->NRr+ir], Beam);
            }
        }
    }
}

/**
 * Bucker's Simple Gaussian Beams in Cartesian coordinates
 */
HOST_DEVICE inline void Step_InfluenceSGB(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray,
    int32_t is, cpx *u,
    const Position *Pos, const BeamStructure *Beam)
{
    real w;
    vec2 x, rayt;
    cpx tau;
    
    const real beta = RC(0.98); // Beam Factor
    real a = RC(-4.0) * STD::log(beta) / SQ(inflray.Dalpha);
    real cn = inflray.Dalpha * STD::sqrt(a / M_PI);
    real rA = point0.x.x;
    int32_t ir = 0;
    
    real rB = point1.x.x;
    
    real q = point0.q.x;
    IncPhaseIfCaustic(inflray, q);
    inflray.qOld = q;
    
    // Loop over bracketted receiver ranges
    // LP: BUG: This does not check for "bracketted" receivers, only that the
    // end of the ray segment is to the right of the receiver. w may be an
    // arbitrarily low negative number, which can lead to garbage extrapolated
    // values for the other variables. See the BUG below.
    while(STD::abs(rB - rA) > RC(1.0e3) * spacing(rA) && rB > Pos->Rr[ir]){
        w    = (Pos->Rr[ir] - rA) / (rB - rA);
        x    = point0.x   + w * (point1.x   - point0.x);
        rayt = point0.t   + w * (point1.t   - point0.t);
        q    = point0.q.x + w * (point1.q.x - point0.q.x);
        tau  = point0.tau + w * (point1.tau - point0.tau);
        
        // following is incorrect because ray doesn't always use a step of deltas
        // LP: The while ignores extremely small steps, but those small steps
        // still increment is, so the later ray segments still treat it as if
        // all steps leading up to them were of size deltas.
        real sint = ((real)is + w) * Beam->deltas;
        
        IncPhaseIfCaustic(inflray, q);
        
        for(int32_t iz=0; i<inflray.NRz_per_range; ++i){
            real deltaz = Pos->Rz[iz] - x.y; // ray to rcvr distance
            //LP: This is commented out, but seems very important to have.
            //Though with all the other bugs, perhaps this is not the main issue.
            //real Adeltaz = STD::abs(deltaz);
            //if(Adeltaz < inflray.RadMax){
                if(Beam->RunType[0] == 'E'){ // eigenrays
                    print("Eigenrays not yet supported\n");
                    bail();
                    //WriteRay2D(inflray.SrcDeclAngle, is);
                }else{ // coherent TL
                    // LP: BUG: It may be incoherent or semi-coherent.
                    real cpa = STD::abs(deltaz * (rB - rA)) / STD::sqrt(SQ(rB - rQ) 
                        + SQ(point1.x.y - point0.x.y));
                    real ds = STD::sqrt(SQ(deltaz) - SQ(cpa));
                    real sx1 = sint + ds;
                    real thet = STD::atan(cpa / sx1);
                    real delay = tau + rayt.y * deltaz;
                    cpx contri = inflray.Ratio1 * cn * point1.Amp * STD::exp(-a * SQ(thet) -
                        J * (inflray.omega * delay - point1.Phase - inflray.phase)) / STD::sqrt(sx1);
                    AtomicAddCpx(&u[iz*Pos->NRr + ir], contri);
                }
            //}
        }
        
        inflray.qOld = q;
        ++ir;
        // LP: BUG: In the FORTRAN this is also return, but it is within the
        // loop over the ray steps. This means that once the ray is to the right
        // of all receivers, it ignores the rest of the ray. This is wrong for
        // two reasons. First, it assumes rays always travel to the right, which
        // is not always true for certain bathymetries. Second, the influence is
        // only computed when the ray is to the right of a receiver. If there
        // are receviers at different range positions, the left ones will
        // receive contributions for every step the ray is to the right. How
        // many of these steps happen depends on where the right receivers are,
        // which means the computed signal at one receiver depends on the
        // position of another receiver, which is absurd.
        if(ir >= Pos->NRr) return;
    }
}

HOST_DEVICE inline void Step_Influence(
    const ray2DPt &point0, const ray2DPt &point1, InfluenceRayInfo &inflray, 
    int32_t is, cpx *u,
    const BdryType *Bdry, const SSPStructure *ssp, int32_t &iSegz, int32_t &iSegr, 
    const Position *Pos, const BeamStructure *Beam)
{
    
    switch(Beam->Type[0]){
    case 'R': Step_InfluenceCervenyRayCen(
            point0, point1, inflray, is, u, Bdry, Pos, Beam
        ); break;
    case 'C': Step_InfluenceCervenyCart(
            point0, point1, inflray, is, u, Bdry, ssp, iSegz, iSegr, Pos, Beam
        ); break;
    case 'g': Step_InfluenceGeoHatRayCen(
            point0, point1, inflray, is, u, Pos, Beam
        ); break;
    case 'S': Step_InfluenceSGB(
            point0, point1, inflray, is, u, Pos, Beam
        ); break;
    case 'B':
    default:  Step_InfluenceGeoHatOrGaussianCart(Beam->Type[0] == 'B',
            point0, point1, inflray, is, u, Pos, Beam
        );
    }
}
