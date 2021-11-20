#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "refcoef.hpp"
#include "ssp.hpp"
#include "attenuation.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"

void setup(std::string FileRoot, 
    std::ofstream &PRTFile, LDOFile &RAYFile, std::ofstream &ARRFile,
    std::string &Title, real &fT,
    BdryType *&Bdry, BdryInfo *&bdinfo, ReflectionInfo *&refl, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, FreqInfo *&freqinfo, 
    BeamStructure *&Beam, BeamInfo *&beaminfo);
    
