/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
Copyright (C) 2021-2022 The Regents of the University of California
c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter

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
#include "common.hpp"
#include "raymode.hpp"
#include "tlmode.hpp"
#include "eigenrays.hpp"
#include "arrivals.hpp"

int main(int argc, char **argv)
{
    std::string FileRoot;
    bool singlethread = false;
    for(int32_t i=1; i<argc; ++i){
        std::string s = argv[i];
        if(argv[i][0] == '-'){
            if(s.length() >= 2 && argv[i][1] == '-'){ //two dashes
                s = s.substr(1);
            }
            if(s == "-1" || s == "-singlethread"){
                singlethread = true;
            }else{
                std::cout << "Unknown command-line option \"" << s << "\"\n";
                std::abort();
            }
        }else{
            if(FileRoot.empty()){
                FileRoot = s;
            }else{
                std::cout << "Intepreting both \"" << FileRoot << "\" and \"" << 
                    s << "\" as FileRoot, error\n";
                std::abort();
            }
        }
    }
    if(FileRoot.empty()){
        std::cout << "Must provide FileRoot as command-line parameter\n";
        std::abort();
    }
    
    bhc::bhcParams params;
    bhc::bhcOutputs outputs;
    bhc::setup(FileRoot.c_str(), nullptr, params, outputs);
    
    bhc::Stopwatch sw;
    sw.tick();
    bhc::run(params, outputs, singlethread);
    sw.tock();
    
    char r = params.Beam->RunType[0];
    if(r == 'R'){
        // Ray mode
        bhc::FinalizeRayMode(outputs.rayinfo, FileRoot, params);
    }else if(r == 'C' || r == 'S' || r == 'I'){
        // TL mode
        bhc::FinalizeTLMode(FileRoot, params, outputs);
    }else if(r == 'E'){
        // Eigenrays mode
        bhc::FinalizeEigenMode(params, outputs, FileRoot, singlethread);
    }else if(r == 'A' || r == 'a'){
        // Arrivals mode
        bhc::FinalizeArrivalsMode(outputs.arrinfo, params.Pos, params.freqinfo,
            params.Beam, FileRoot, false);
    }else{
        std::cout << "Invalid RunType " << r << "\n";
        std::abort();
    }
    
    bhc::finalize(params, outputs);
}
