#pragma once
#include "common.hpp"

void ReadEnvironment(const std::string &FileRoot, std::ostream &PRTFile,
    std::string &Title, real &fT, BdryType *Bdry, SSPStructure *ssp, AttenInfo *atten, 
    Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, BeamStructure *Beam,
    HSInfo &RecycledHS);
