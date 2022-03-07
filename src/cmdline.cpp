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
    OpenPRTFile(FileRoot, PRTFile);
    
    bhcParams params;
    bhcOutputs outputs;
    setup(FileRoot, PRTFile, params, outputs);
    
    Stopwatch sw;
    sw.tick();
    run(PRTFile, params, outputs, singlethread);
    sw.tock();
    
    char r = params.Beam->RunType[0];
    if(r == 'R'){
        // Ray mode
        FinalizeRayMode(outputs.rayinfo, FileRoot, params);
    }else if(r == 'C' || r == 'S' || r == 'I'){
        // TL mode
        FinalizeTLMode(FileRoot, params, outputs);
    }else if(r == 'E'){
        // Eigenrays mode
        FinalizeEigenMode(params, outputs, FileRoot, singlethread);
    }else if(r == 'A' || r == 'a'){
        // Arrivals mode
        FinalizeArrivalsMode(outputs.arrinfo, params.Pos, params.freqinfo,
            params.Beam, FileRoot, false);
    }else{
        std::cout << "Invalid RunType " << r << "\n";
        std::abort();
    }
}
