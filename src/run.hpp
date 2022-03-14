#pragma once
#include "common.hpp"
#include "raymode.hpp"
#include "tlmode.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"

namespace bhc {

HOST_DEVICE inline int32_t GetNumJobs(const Position *Pos, const AnglesStructure *Angles)
{
    return Pos->NSz * (Angles->iSingle_alpha >= 1 ? 1 : Angles->Nalpha);
}

/**
 * Returns whether the job should continue.
 * `is` changed to `isrc` because `is` is used for steps
 */
HOST_DEVICE inline bool GetJobIndices(int32_t &isrc, int32_t &ialpha, int32_t job,
    const Position *Pos, const AnglesStructure *Angles)
{
    if(Angles->iSingle_alpha >= 1){
        isrc = job;
        ialpha = Angles->iSingle_alpha - 1; //iSingle_alpha is 1-indexed because how defined in env file
    }else{
        isrc = job / Angles->Nalpha;
        ialpha = job % Angles->Nalpha;
    }
    return (isrc < Pos->NSz);
}

inline void InitSelectedMode(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread)
{
    // Common
    // irregular or rectilinear grid
    params.Pos->NRz_per_range = (params.Beam->RunType[4] == 'I') ? 1 : params.Pos->NRz;
    
    // Mode specific
    if(params.Beam->RunType[0] == 'R'){
        InitRayMode(outputs.rayinfo, params);
    }else if(params.Beam->RunType[0] == 'C' || params.Beam->RunType[0] == 'S' || params.Beam->RunType[0] == 'I'){
        InitTLMode(outputs.uAllSources, params.Pos);
    }else if(params.Beam->RunType[0] == 'E'){
        InitEigenMode(outputs.eigen);
    }else if(params.Beam->RunType[0] == 'A' || params.Beam->RunType[0] == 'a'){
        InitArrivalsMode(outputs.arrinfo, singlethread, params.Pos, PRTFile);
    }else{
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        std::abort();
    }
}

void run_cxx(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread);

}
