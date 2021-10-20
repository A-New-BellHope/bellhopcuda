#pragma once
#include "common.hpp"

HOST_DEVICE inline void SubTab(real *x, int32_t Nx)
{
    if(Nx >= 3){
        if(x[2] == RC(-999.9)){ // testing for equality here is dangerous
            if(x[1] == RC(-999.9)) x[1] = x[0];
            real deltax = (x[1] - x[0]) / (real)(Nx - 1);
            real x0 = x[0];
            for(int32_t i=0; i<Nx; ++i) x[i] = x0 + (real)i * deltax;
        }
    }
}
