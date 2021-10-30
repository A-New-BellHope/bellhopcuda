#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "attenuation.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

void ReadEnvironment(const std::string &FileRoot, std::ofstream &PRTFile,
    std::string &Title, real &fT, BdryType *Bdry, SSPStructure *ssp, AttenInfo *atten, 
    Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, BeamStructure *Beam);
