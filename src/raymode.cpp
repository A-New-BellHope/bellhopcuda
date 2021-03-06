/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

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
#include "raymode.hpp"
#include "run.hpp"

namespace bhc {

void OpenRAYFile(LDOFile &RAYFile, std::string FileRoot, bool ThreeD, 
    const bhcParams &params)
{
    switch(params.Beam->RunType[0]){
    case 'R':
    case 'E':
        // Ray trace or Eigenrays
        break;
    default:
        std::cout << "OpenRAYFile not in ray trace or eigenrays mode\n";
        std::abort();
    }
    RAYFile.open(FileRoot + ".ray");
    RAYFile << params.Title << '\n';
    RAYFile << params.freqinfo->freq0 << '\n';
    RAYFile << params.Pos->NSx << params.Pos->NSy << params.Pos->NSz << '\n';
    RAYFile << params.Angles->Nalpha << params.Angles->Nbeta << '\n';
    RAYFile << params.Bdry->Top.hs.Depth << '\n';
    RAYFile << params.Bdry->Bot.hs.Depth << '\n';
    RAYFile << (ThreeD ? "xyz" : "rz") << '\n';
}


/**
 * Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
 * Write to RAYFile.
 * 
 * During an eigenray calculation, subsets of the full ray may be passed
 * These have lengths Nsteps1 vs. Nsteps for the entire ray
 * 
 * The 2D version is for ray traces in (r,z) coordinates
 * 
 * alpha0: take-off angle of this ray
 */
void WriteRay2D(real alpha0, int32_t Nsteps1, LDOFile &RAYFile,
    const BdryType *Bdry, ray2DPt *ray2D)
{
    // compression
    
    constexpr int32_t MaxNRayPoints = 500000; // this is the maximum length of the ray vector that is written out
    int32_t n2 = 1;
    int32_t iSkip = bhc::max(Nsteps1 / MaxNRayPoints, 1);
    
    for(int32_t is=1; is<Nsteps1; ++is){
        // ensure that we always write ray points near bdry reflections (works only for flat bdry)
        if(bhc::min(Bdry->Bot.hs.Depth - ray2D[is].x.y, ray2D[is].x.y - Bdry->Top.hs.Depth) < FL(0.2) ||
                (is % iSkip) == 0 || is == Nsteps1-1){
            ++n2;
            ray2D[n2-1].x = ray2D[is].x;
        }
    }
    
    // write to ray file
    
    RAYFile << alpha0 << '\n';
    RAYFile << n2 << ray2D[Nsteps1-1].NumTopBnc << ray2D[Nsteps1-1].NumBotBnc << '\n';
    
    for(int32_t is=0; is<n2; ++is){
        RAYFile << ray2D[is].x << '\n';
    }
}

void InitRayMode(RayInfo *rayinfo, const bhcParams &params)
{
    rayinfo->NRays = GetNumJobs(params.Pos, params.Angles);
    rayinfo->MaxPoints = bhc::min((uint32_t)MaxN * (uint32_t)rayinfo->NRays, 100000000u);
    rayinfo->NPoints = 0;
    checkallocate(rayinfo->raymem, rayinfo->MaxPoints);
    checkallocate(rayinfo->results, rayinfo->NRays);
    memset(rayinfo->results, 0, rayinfo->NRays * sizeof(RayResult)); // Clear because will check pointers
}

void FinalizeRayMode(RayInfo *rayinfo, std::string FileRoot, const bhcParams &params)
{
    LDOFile RAYFile;
    OpenRAYFile(RAYFile, FileRoot, false, params);
    for(int r=0; r<rayinfo->NRays; ++r){
        RayResult *res = &rayinfo->results[r];
        if(res->ray2D == nullptr) continue;
        WriteRay2D(res->SrcDeclAngle, res->Nsteps, RAYFile, params.Bdry, res->ray2D);
    }
}

}
