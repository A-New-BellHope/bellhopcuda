/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/

#include <iostream>

// This define must be set before including the header if you're using the DLL
// version on Windows, and it must NOT be set if you're using the static library
// version on Windows. If you're not on Windows, it doesn't matter either way.
#define BHC_DLL_IMPORT 1
#include <bhc/bhc.hpp>

/*
std::ostream& operator<<(std::ostream& out, const bhc::rayPt<false>& x) {
    out << x.NumTopBnc << " "
        << x.NumBotBnc << " "
        << x.x.x << " " << x.x.y << " "
        << x.t.x << " " << x.t.y << " "
        << x.p.x << " " << x.p.y << " "
        << x.q.x << " " << x.q.y << " "
        << x.c << " "
        << x.Amp << " " << x.Phase << " "
        << x.tau;
    return out;
}
*/
std::ostream &operator<<(std::ostream &out, const bhc::cpxf &x)
{
    out << "(" << x.real() << ", " << x.imag() << ")";
    return out;
}

void OutputCallback(const char *message)
{
    std::cout << "Out: " << message << std::endl << std::flush;
}

void PrtCallback(const char *message) { std::cout << message << std::flush; }

int main()
{
    bhc::bhcParams<false> params;
    bhc::bhcOutputs<false, false> outputs;
    bhc::bhcInit init;
    init.FileRoot       = nullptr;
    init.outputCallback = OutputCallback;
    init.prtCallback    = PrtCallback;

    bhc::setup(init, params, outputs);
    bhc::echo(params);
    bhc::run(params, outputs);

    std::cout << "Field:\n";
    for(int iz = 0; iz < params.Pos->NRz; ++iz) {
        for(int ir = 0; ir < params.Pos->NRr; ++ir) {
            std::cout << outputs.uAllSources[iz * params.Pos->NRr + ir] << " ";
        }
        std::cout << "\n";
    }
    /*
    std::cout << "\n" << outputs.rayinfo->NRays << " rays:\n";
    for(int r=0; r<outputs.rayinfo->NRays; ++r){
        std::cout << "\nRay " << r
            << ", " << outputs.rayinfo->results[r].Nsteps
            << "steps, SrcDeclAngle = " << outputs.rayinfo->results[r].SrcDeclAngle
            << ":\n";
        for(int s=0; s<outputs.rayinfo->results[r].Nsteps; ++s){
            std::cout << outputs.rayinfo->results[r].ray[s] << "\n";
        }
    }
    */

    bhc::finalize(params, outputs);
    return 0;
}
