#pragma once
#include "common.hpp"

struct Arrival {
    int32_t NTopBnc, NBotBnc;
    float SrcDeclAngle, SrcAzimAngle, RcvrDeclAngle, RcvrAzimAngle, a, Phase;
    cpxf delay;
};

/**
 * Variables for arrival information
 */
struct ArrInfo {
    Arrival *Arr;
    int32_t *NArr;
    int32_t MaxNArr;
};

constexpr int32_t ArrivalsStorage = 20000000, MinNArr = 10;

/**
 * ADDs the amplitude and delay for an ARRival into a matrix of same.
 * Extra logic included to keep only the strongest arrivals.
 */
HOST_DEVICE inline void AddArr(real omega, int32_t isrc, int32_t id, int32_t ir,
    real Amp, real Phase, cpx delay, real SrcDeclAngle, real RcvrDeclAngle,
    int32_t NumTopBnc, int32_t NumBotBnc,
    const ArrInfo *arrinfo, const Position *Pos)
{
    const float PhaseTol = FL(0.05); // arrivals with essentially the same phase are grouped into one
    
    int32_t base = (isrc * Pos->NRz_per_range + id) * Pos->NRr + ir;
    Arrival *baseArr = &arrinfo->Arr[base * arrinfo->MaxNArr];
    int32_t *baseNArr = &arrinfo->NArr[base];
    
    int32_t Nt = *baseNArr; // # of arrivals
    bool NewRay = true;
    
    // Is this the second bracketting ray of a pair?
    // If so, we want to combine the arrivals to conserve space.
    // (test this by seeing if the arrival time is close to the previous one)
    // (also need that the phase is about the same to make sure surface and direct paths are not joined)
    
    if(Nt >= 1){
        if(omega * STD::abs(delay - baseArr[Nt-1].delay) < PhaseTol &&
            STD::abs(baseArr[Nt-1].phase - Phase) < PhaseTol) NewRay = false;
    }
    
    if(NewRay){
        int32_t iArr;
        if(Nt >= arrinfo->MaxNArr){ // space [LP: NOT] available to add an arrival?
            iArr = TODO; // no: replace weakest arrival
            if(Amp <= baseArr[iArr].a) return;
        }else{
            iArr = Nt;
            *baseNArr = Nt + 1; // # of arrivals
        }
        baseArr[iArr].A             = (float)Amp; // amplitude
        baseArr[iArr].Phase         = (float)Phase; // phase
        baseArr[iArr].delay         = Cpx2Cpxf(delay); // delay time
        baseArr[iArr].SrcDeclAngle  = (float)SrcDeclAngle; // angle
        baseArr[iArr].RcvrDeclAngle = (float)RcvrDeclAngle; // angle [LP: :( ]
        baseArr[iArr].NTopBnc       = NumTopBnc; // Number of top    bounces
        baseArr[iArr].NBotBnc       = NumBotBnc; //   "       bottom
    }else{ // not a new ray
        // PhaseArr[<base> + Nt-1] = PhaseArr[<base> + Nt-1] // LP: ???
        
        // calculate weightings of old ray information vs. new, based on amplitude of the arrival
        float AmpTot = baseArr[Nt-1].a + (float)Amp;
        float w1 = baseArr[Nt-1].a / AmpTot;
        float w2 = (float)Amp / AmpTot;
        
        baseArr[Nt-1].delay         = w1 * baseArr[Nt-1].delay         + w2 * Cpx2Cpxf(delay); // weighted sum
        baseArr[Nt-1].a             = AmpTot;
        baseArr[Nt-1].SrcDeclAngle  = w1 * baseArr[Nt-1].SrcDeclAngle  + w2 * (float)SrcDeclAngle;
        baseArr[Nt-1].RcvrDeclAngle = w1 * baseArr[Nt-1].RcvrDeclAngle + w2 * (float)RcvrDeclAngle;
    }
}

inline void InitArrivalsMode(ArrInfo *arrinfo, 
    const Position *Pos, const BeamStructure *Beam, LDOFile &PRTFile)
{
    int32_t nzr = Pos->NRz_per_range * Pos->NRr;
    arrinfo->MaxNArr = math::max(ArrivalsStorage / nzr, MinNArr);
    PRTFile << "\n( Maximum # of arrivals = " << arrinfo->MaxNArr << ")\n";
    arrinfo->Arr = allocate<Arrival>(Pos->NSz * nzr * arrinfo->MaxNArr);
    arrinfo->NArr = allocate<int32_t>(Pos->NSz * nzr);
    memset(arrinfo->NArr, 0, Pos->NSz * nzr * sizeof(int32_t));
    memset(arrinfo->Narr, 0, Pos->NSz * nzr * arrinfo->MaxNArr * sizeof(Arrival));
}

void WriteArrivals(const ArrInfo *arrinfo, const Position *Pos,
    const BeamStructure *Beam, std::string FileRoot, bool ThreeD);
