#pragma once
#include "common.hpp"

struct AnglesStructure {
    int32_t Nalpha = 0, Nbeta = 1, iSingle_alpha = -1, iSingle_beta = -1;
    real Dalpha, Dbeta;
    real *alpha, *beta;
};
