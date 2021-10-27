#pragma once
#include "common.hpp"
#include "sourcereceiver.hpp"
#include "boundary.hpp"
#include "ssp.hpp"
#include "attenuation.hpp"
#include "angles.hpp"
#include "beampattern.hpp"
#include "refcoef.hpp"

void setup(int argc, char **argv, std::string &Title, real &fT,
    FreqInfo *&freqinfo, BdryType *&Bdry, BdryInfo *&binfo, SSPStructure *&ssp,
    AttenInfo *&atten, Position *&Pos, AnglesStructure *&Angles, BeamStructure *&Beam,
    BeamInfo *&beaminfo, ReflectionInfo *&refl);
