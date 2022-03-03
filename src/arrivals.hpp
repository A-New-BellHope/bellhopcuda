#pragma once
#include "common.hpp"

/**
 * Is this the second step of a pair (on the same ray)?
 * If so, we want to combine the arrivals to conserve space.
 * (test this by seeing if the arrival time is close to the previous one)
 * (also need that the phase is about the same to make sure surface and direct paths are not joined)
 */
HOST_DEVICE inline bool IsSecondStepOfPair(real omega, real Phase, cpx delay,
    Arrival *baseArr, int32_t Nt)
{
    const float PhaseTol = FL(0.05); // arrivals with essentially the same phase are grouped into one
    return Nt >= 1 && 
        omega * STD::abs(delay - Cpxf2Cpx(baseArr[Nt-1].delay)) < PhaseTol &&
        STD::abs(baseArr[Nt-1].Phase - Phase) < PhaseTol;
}

/**
 * ADDs the amplitude and delay for an ARRival into a matrix of same.
 * Extra logic included to keep only the strongest arrivals.
 */
HOST_DEVICE inline void AddArr(real omega, int32_t isrc, int32_t id, int32_t ir,
    real Amp, real Phase, cpx delay, real SrcDeclAngle, real RcvrDeclAngle,
    int32_t NumTopBnc, int32_t NumBotBnc,
    const ArrInfo *arrinfo, const Position *Pos)
{
    size_t base = (isrc * Pos->NRz_per_range + id) * Pos->NRr + ir;
    Arrival *baseArr = &arrinfo->Arr[base * arrinfo->MaxNArr];
    int32_t *baseNArr = &arrinfo->NArr[base];
    int32_t Nt;
    
    if(arrinfo->singlethread){
        // LP: Besides this algorithm being virtually impossible to adapt for
        // multithreading, there is a more fundamental BUG with it: it does
        // not have any handling for the case where the current arrival is the
        // second step of a pair, but the storage is full and the first step of
        // the pair got stored in some slot other than the last one. It will
        // only consider the last ray in storage as a candidate for the first
        // half of the pair. This means that whether a given pair is
        // successfully paired or not depends on the number of arrivals before
        // that pair, which depends on the order rays were computed, which is
        // arbitrary and non-physical.
    
        Nt = *baseNArr; // # of arrivals
        
        if(!IsSecondStepOfPair(omega, Phase, delay, baseArr, Nt)){
            int32_t iArr;
            if(Nt >= arrinfo->MaxNArr){ // space not available to add an arrival?
                // replace weakest arrival
                iArr = -1;
                real weakest = Amp;
                for(Nt=0; Nt<arrinfo->MaxNArr; ++Nt){
                    if(baseArr[Nt].a < weakest){
                        weakest = baseArr[Nt].a;
                        iArr = Nt;
                    }
                }
                if(iArr < 0) return; // LP: current arrival is weaker than all stored
            }else{
                iArr = Nt;
                *baseNArr = Nt + 1; // # of arrivals
            }
            baseArr[iArr].a             = (float)Amp; // amplitude
            baseArr[iArr].Phase         = (float)Phase; // phase
            baseArr[iArr].delay         = Cpx2Cpxf(delay); // delay time
            baseArr[iArr].SrcDeclAngle  = (float)SrcDeclAngle; // launch angle from source
            baseArr[iArr].RcvrDeclAngle = (float)RcvrDeclAngle; // angle ray reaches receiver
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
        
    }else{
        // LP: For multithreading mode, some mutex scheme would be needed to
        // guarantee correct access to previously written data, which would
        // destroy the performance on GPU. So just write the first 
        // arrinfo->MaxNArr arrivals and give up.
        Nt = AtomicFetchAdd(baseNArr, 1);
        if(Nt >= arrinfo->MaxNArr) return;
        baseArr[Nt].a             = (float)Amp; // amplitude
        baseArr[Nt].Phase         = (float)Phase; // phase
        baseArr[Nt].delay         = Cpx2Cpxf(delay); // delay time
        baseArr[Nt].SrcDeclAngle  = (float)SrcDeclAngle; // launch angle from source
        baseArr[Nt].RcvrDeclAngle = (float)RcvrDeclAngle; // angle ray reaches receiver
        baseArr[Nt].NTopBnc       = NumTopBnc; // Number of top    bounces
        baseArr[Nt].NBotBnc       = NumBotBnc; //   "       bottom
    }
}

inline void InitArrivalsMode(ArrInfo *arrinfo, bool singlethread,
    const Position *Pos, const BeamStructure *Beam, std::ofstream &PRTFile)
{
    arrinfo->singlethread = singlethread;
    size_t nzr = Pos->NRz_per_range * Pos->NRr;
    const size_t ArrivalsStorage = singlethread ? 20000000 : 50000000;
    const size_t MinNArr = 10;
    arrinfo->MaxNArr = math::max(ArrivalsStorage / nzr, MinNArr);
    PRTFile << "\n( Maximum # of arrivals = " << arrinfo->MaxNArr << ")\n";
    arrinfo->Arr  = allocate<Arrival>((size_t)Pos->NSz * nzr * (size_t)arrinfo->MaxNArr);
    arrinfo->NArr = allocate<int32_t>((size_t)Pos->NSz * nzr);
    memset(arrinfo->Arr,  0, (size_t)Pos->NSz * nzr * (size_t)arrinfo->MaxNArr * sizeof(Arrival));
    memset(arrinfo->NArr, 0, (size_t)Pos->NSz * nzr * sizeof(int32_t));
}

void FinalizeArrivalsMode(const ArrInfo *arrinfo, const Position *Pos,
    const FreqInfo *freqinfo, const BeamStructure *Beam, std::string FileRoot, bool ThreeD);
