#pragma once
#include "common.hpp"

HOST_DEVICE inline void RecordEigenHit(int32_t ir, int32_t iz, 
    int32_t isrc, int32_t ialpha, int32_t is, EigenInfo *eigen)
{
    uint32_t mi = AtomicFetchAdd(&eigen->neigen, 1u);
    if(mi >= eigen->memsize) return;
    // printf("Eigenray hit %d ir %d iz %d isrc %d ialpha %d is %d\n",
    //     mi, ir, iz, isrc, ialpha, is);
    eigen->hits[mi].ir = ir;
    eigen->hits[mi].iz = iz;
    eigen->hits[mi].isrc = isrc;
    eigen->hits[mi].ialpha = ialpha;
    eigen->hits[mi].is = is;
}

inline void InitEigenMode(EigenInfo *eigen)
{
    static const uint32_t maxhits = 1000000u;
    eigen->neigen = 0;
    eigen->memsize = maxhits;
    eigen->hits = allocate<EigenHit>(maxhits);
}

void FinalizeEigenMode(
    const BdryType *Bdry, const BdryInfo *bdinfo, const ReflectionInfo *refl,
    const SSPStructure *ssp, const Position *Pos, const AnglesStructure *Angles,
    const FreqInfo *freqinfo, const BeamStructure *Beam, const BeamInfo *beaminfo,
    const EigenInfo *eigen, LDOFile &RAYFile, bool singlethread);
