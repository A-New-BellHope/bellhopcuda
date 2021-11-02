#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "refcoef.hpp"
#include "ssp.hpp"
#include "attenuation.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

void setup(int argc, char **argv, 
    std::ofstream &PRTFile, std::ofstream &RAYFile, std::ofstream &ARRFile,
    std::string &Title, real &fT,
    BdryType *&Bdry, BdryInfo *&bdinfo, ReflectionInfo *&refl, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, FreqInfo *&freqinfo, 
    BeamStructure *&Beam, BeamInfo *&beaminfo);
    
void core_setup(std::ofstream &PRTFile, const real &fT,
    const BdryType *Bdry, const BdryInfo *bdinfo, const AttenInfo *atten, 
    AnglesStructure *Angles, const FreqInfo *freqinfo, BeamStructure *Beam/*, 
    InfluenceInfo *inflinfo, ArrivalsInfo *arrinfo*/);
