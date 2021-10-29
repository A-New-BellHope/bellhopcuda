#include "common.hpp"

/**
 * Write appropriate header information
 */
void OpenOutputFiles(std::string FileRoot, bool ThreeD, std::string Title,
    BdryType *Bdry, Position *Pos, AnglesStructure *Angles, FreqInfo *freqinfo, 
    BeamStructure *Beam,
    std::ofstream &RAYFile, std::ofstream &ARRFile)
{
    real atten;
    std::string PlotType;
    
    switch(Beam->RunType[0]){
    case 'R':
    case 'E':
        // Ray trace or Eigenrays
        RAYFile.open(WithExtension(FileRoot, ".ray"));
        RAYFile << "'" << Title << "'\n";
        RAYFile << freqinfo->freq0 << "\n";
        RAYFile << Pos->NSx << " " << Pos->NSy << " " << Pos->NSz << "\n";
        RAYFile << Angles->Nalpha << " " << Angles->Nbeta << "\n";
        RAYFile << Bdry->Top.hs.Depth << "\n";
        RAYFile << Bdry->Bot.hs.Depth << "\n";
        
        RAYFile << (ThreeD ? "'xyz'" : "'rz'") << "\n";
        break;
    case 'A':
        // arrival file in ascii format
        ARRFile.open(WithExtension(FileRoot, ".arr"));
        
        ARRFile << (ThreeD ? "'3D'" << "'2D'") << "\n";
        
        ARRFile << freqinfo->freq0 << "\n";
        
        // write source locations
        
        if(ThreeD){
            ARRFile << Pos->NSx; for(int32_t i=0; i<Pos->NSx; ++i) ARRFile << " " << Pos->Sx[i]; ARRFile << "\n";
            ARRFile << Pos->NSy; for(int32_t i=0; i<Pos->NSy; ++i) ARRFile << " " << Pos->Sy[i]; ARRFile << "\n";
        }
        ARRFile << Pos->NSz; for(int32_t i=0; i<Pos->NSz; ++i) ARRFile << " " << Pos->Sz[i]; ARRFile << "\n";
        
        // write receiver locations
        ARRFile << Pos->NRz; for(int32_t i=0; i<Pos->NRz; ++i) ARRFile << " " << Pos->Rz[i]; ARRFile << "\n";
        ARRFile << Pos->NRr; for(int32_t i=0; i<Pos->NRr; ++i) ARRFile << " " << Pos->Rr[i]; ARRFile << "\n";
        if(ThreeD){
            ARRFile << Pos->Ntheta; for(int32_t i=0; i<Pos->Ntheta; ++i) ARRFile << " " << Pos->theta[i]; ARRFile << "\n";
        }
        break;
    case 'a':
        // arrival file in binary format
        ARRFile.open(WithExtension(FileRoot, ".arr"));
        
        ARRFile.write((ThreeD ? "'3D'" << "'2D'"), 4);
        
        float freq = freqinfo->freq0; ARRFile.write(&freq, 4);
        
        // write source locations
        
        #define RS sizeof(real)
        if(ThreeD){
            ARRFile.write(&Pos->NSx, 4); for(int32_t i=0; i<Pos->NSx; ++i) ARRFile.write(&Pos->Sx[i], RS);
            ARRFile.write(&Pos->NSy, 4); for(int32_t i=0; i<Pos->NSy; ++i) ARRFile.write(&Pos->Sy[i], RS);
        }
        ARRFile.write(&Pos->NSz, 4); for(int32_t i=0; i<Pos->NSz; ++i) ARRFile.write(&Pos->Sz[i], RS);
        
        // write receiver locations
        ARRFile.write(&Pos->NRz, 4); for(int32_t i=0; i<Pos->NRz; ++i) ARRFile.write(&Pos->Rz[i], RS);
        ARRFile.write(&Pos->NRr, 4); for(int32_t i=0; i<Pos->NRr; ++i) ARRFile.write(&Pos->Rr[i], RS);
        if(ThreeD){
            ARRFile.write(&Pos->Ntheta, 4); for(int32_t i=0; i<Pos->Ntheta; ++i) ARRFile.write(&Pos->theta[i], RS);
        }
        #undef RS
        break;
    default:
        atten = RC(0.0);
        
        // following to set PlotType has alread been done in READIN if that was used for input
        switch(Beam->RunType[4]){
        case 'R':
            PlotType = "rectilin  "; break;
        case 'I':
            PlotType = "irregular "; break;
        default:
            PlotType = "rectilin  ";
        }
        WriteHeader(WithExtension(FileRoot, ".shd"), Title, freqinfo->freq0, atten, PlotType);
    }
}
