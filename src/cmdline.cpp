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
    
    std::ofstream PRTFile;
    OpenPRTFile(PRTFile);
    
    bhcParams params;
    setup(FileRoot, PRTFile, params);
    
    Stopwatch sw;
    sw.tick();
    bchOutputs outputs;
    run(PRTFile, params, outputs, singlethread);
    sw.tock();
    
    if(Beam->RunType[0] == 'R'){
        // Ray mode
        FinalizeRayMode(outputs.rayinfo, FileRoot, params);
    }else if(Beam->RunType[0] == 'C' || Beam->RunType[0] == 'S' || Beam->RunType[0] == 'I'){
        // TL mode
        FinalizeTLMode(outputs.uAllSources, SHDFile, ssp, Pos, Angles, freqinfo, Beam);
    }else if(Beam->RunType[0] == 'E'){
        // Eigenrays mode
        FinalizeEigenMode(params, outputs, FileRoot, singlethread);
    }else if(Beam->RunType[0] == 'A' || Beam->RunType[0] == 'a'){
        // Arrivals mode
        FinalizeArrivalsMode(outputs.arrinfo, Pos, freqinfo, Beam, FileRoot, false);
    }else{
        std::cout << "Invalid RunType " << Beam->RunType[0] << "\n";
        std::abort();
    }
}
