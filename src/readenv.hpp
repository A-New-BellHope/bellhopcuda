#pragma once
#include "common.hpp"

namespace bhc {

void ReadEnvironment(const std::string &FileRoot, std::ostream &PRTFile,
    char (&Title)[80], real &fT, BdryType *Bdry, SSPStructure *ssp, AttenInfo *atten, 
    Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, BeamStructure *Beam,
    HSInfo &RecycledHS);

}
