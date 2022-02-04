#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "refcoef.hpp"
#include "ssp.hpp"
#include "attenuation.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"
#include "eigenrays.hpp"

void setup(std::string FileRoot, 
    std::ofstream &PRTFile, LDOFile &RAYFile, DirectOFile &SHDFile,
    std::string &Title, real &fT,
    BdryType *&Bdry, BdryInfo *&bdinfo, ReflectionInfo *&refl, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, FreqInfo *&freqinfo, 
    BeamStructure *&Beam, BeamInfo *&beaminfo, EigenInfo *&eigen);
    
/**
 * Call after changing parameters (or setup), and before the init for your
 * desired mode.
 */
inline void InitCommon(Position *Pos, const BeamStructure *Beam)
{
    Pos->NRz_per_range = (Beam->RunType[4] == 'I') ? 1 : Pos->NRz; // irregular or rectilinear grid
}
