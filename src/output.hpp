#pragma once
#include "common.hpp"
#include "boundary.hpp"
#include "sourcereceiver.hpp"
#include "angles.hpp"
#include "beams.hpp"
#include "step.hpp"

void WriteRay2D(real alpha0, int32_t Nsteps1, std::ofstream &RAYFile,
    const BdryType *Bdry, ray2DPt *ray2D);

void OpenOutputFiles(std::string FileRoot, bool ThreeD, std::string Title,
    const BdryType *Bdry, const Position *Pos, const AnglesStructure *Angles, 
    const FreqInfo *freqinfo, const BeamStructure *Beam,
    std::ofstream &RAYFile, std::ofstream &ARRFile);
