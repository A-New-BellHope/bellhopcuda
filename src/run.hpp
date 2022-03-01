#pragma once
#include "common.hpp"

void run_cxx(std::ostream &PRTFile, const bhcParams &params, bhcOutputs &outputs,
    bool singlethread);
    
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
        InitTLMode(outputs.uAllSources, params.Pos, params.Beam);
    }else if(params.Beam->RunType[0] == 'E'){
        InitEigenMode(outputs.eigen);
    }else if(params.Beam->RunType[0] == 'A' || params.Beam->RunType[0] == 'a'){
        InitArrivalsMode(outputs.arrinfo, singlethread, params.Pos, params.Beam, PRTFile);
    }else{
        std::cout << "Invalid RunType " << params.Beam->RunType[0] << "\n";
        std::abort();
    }
}
