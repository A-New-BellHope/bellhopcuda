#include "output.hpp"

constexpr int32_t MaxNRayPoints = 500000; // this is the maximum length of the ray vector that is written out

/**
 * Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
 * Write to RAYFile.
 * 
 * During an eigenray calculation, subsets of the full ray may be passed
 * These have lengths Nsteps1 vs. Nsteps for the entire ray
 * 
 * The 2D version is for ray traces in (r,z) coordinates
 * 
 * alpha0: take-off angle of this ray
 */
void WriteRay2D(real alpha0, int32_t Nsteps1, LDOFile &RAYFile,
    const BdryType *Bdry, ray2DPt *ray2D)
{
    // compression
    
    int32_t n2 = 1;
    int32_t iSkip = std::max(Nsteps1 / MaxNRayPoints, 1);
    
    for(int32_t is=1; is<Nsteps1; ++is){
        // ensure that we always write ray points near bdry reflections (works only for flat bdry)
        if(std::min(Bdry->Bot.hs.Depth - ray2D[is].x.y, ray2D[is].x.y - Bdry->Top.hs.Depth) < 0.2 ||
                (is % iSkip) == 0 || is == Nsteps1-1){
            ++n2;
            ray2D[n2-1].x = ray2D[is].x;
        }
    }
    
    // write to ray file
    
    RAYFile << alpha0 << '\n';
    RAYFile << n2 << ray2D[Nsteps1-1].NumTopBnc << ray2D[Nsteps1-1].NumBotBnc << '\n';
    
    for(int32_t is=0; is<n2; ++is){
        RAYFile << ray2D[is].x << '\n';
    }
}

/**
 * Write header to disk file
 * LP: of a SHDFile ("binary `shade' file (SHDFIL) that contains calculated
 * pressure fields")
 *
 * FileName: Name of the file (could be a shade file or a Green's function file)
 * Title: Arbitrary title
 * freq0: Nominal frequency [LP: now in freqinfo]
 * Atten: stabilizing attenuation (for wavenumber integration only)
 * PlotType: [LP: argument description comment present but blank; also never
 * set to "TL" in BELLHOP]
 */
void WriteHeader(DirectOFile &SHDFile, const std::string &FileName, 
    const std::string &Title, float Atten, const std::string &PlotType, 
    const Position *Pos, const FreqInfo *freqinfo)
{
    bool isTL = (PlotType[0] == 'T' && PlotType[1] == 'L');
    
    int32_t LRecl = 84; //4 for LRecl, 80 for Title
    LRecl = std::max(LRecl, 2 * freqinfo->Nfreq * (int32_t)sizeof(freqinfo->freqVec[0]));
    LRecl = std::max(LRecl, Pos->Ntheta * (int32_t)sizeof(Pos->theta[0]));
    if(!isTL){
        LRecl = std::max(LRecl, Pos->NSx * (int32_t)sizeof(Pos->Sx[0]));
        LRecl = std::max(LRecl, Pos->NSy * (int32_t)sizeof(Pos->Sy[0]));
    }
    LRecl = std::max(LRecl, Pos->NSz * (int32_t)sizeof(Pos->Sz[0]));
    LRecl = std::max(LRecl, Pos->NRz * (int32_t)sizeof(Pos->Rz[0]));
    LRecl = std::max(LRecl, Pos->NRr * (int32_t)sizeof(cpxf));
    
    SHDFile.open(FileName, LRecl);
    if(!SHDFile.good()){
        std::cout << "Could not open SHDFile: " << FileName << "\n";
        std::abort();
    }
    LRecl /= 4;
    SHDFile.rec(0); DOFWRITE(SHDFile, &LRecl, 4); DOFWRITE(SHDFile, Title, 80);
    SHDFile.rec(1); DOFWRITE(SHDFile, PlotType, 10);
    SHDFile.rec(2);
        DOFWRITEV(SHDFile, freqinfo->Nfreq);
        DOFWRITEV(SHDFile, Pos->Ntheta);
        DOFWRITEV(SHDFile, Pos->NSx);
        DOFWRITEV(SHDFile, Pos->NSy);
        DOFWRITEV(SHDFile, Pos->NSz);
        DOFWRITEV(SHDFile, Pos->NRz);
        DOFWRITEV(SHDFile, Pos->NRr);
        DOFWRITEV(SHDFile, (float)freqinfo->freq0);
        DOFWRITEV(SHDFile, Atten);
    SHDFile.rec(3); DOFWRITE(SHDFile, freqinfo->freqVec, freqinfo->Nfreq * sizeof(freqinfo->freqVec[0]));
    SHDFile.rec(4); DOFWRITE(SHDFile, Pos->theta, Pos->Ntheta * sizeof(Pos->theta[0]));
    
    if(!isTL){
        SHDFile.rec(5); DOFWRITE(SHDFile, Pos->Sx, Pos->NSx * sizeof(Pos->Sx[0]));
        SHDFile.rec(6); DOFWRITE(SHDFile, Pos->Sy, Pos->NSy * sizeof(Pos->Sy[0]));
    }else{
        SHDFile.rec(5); DOFWRITEV(SHDFile, Pos->Sx[0]); DOFWRITEV(SHDFile, Pos->Sx[Pos->NSx-1]);
        SHDFile.rec(6); DOFWRITEV(SHDFile, Pos->Sy[0]); DOFWRITEV(SHDFile, Pos->Sy[Pos->NSy-1]);
    }
    SHDFile.rec(7); DOFWRITE(SHDFile, Pos->Sz, Pos->NSz * sizeof(Pos->Sz[0]));
    
    SHDFile.rec(8); DOFWRITE(SHDFile, Pos->Rz, Pos->NRz * sizeof(Pos->Rz[0]));
    SHDFile.rec(9); DOFWRITE(SHDFile, Pos->Rr, Pos->NRr * sizeof(Pos->Rr[0]));
}

/**
 * Write appropriate header information
 */
void OpenOutputFiles(std::string FileRoot, bool ThreeD, std::string Title,
    const BdryType *Bdry, const Position *Pos, const AnglesStructure *Angles, 
    const FreqInfo *freqinfo, const BeamStructure *Beam,
    LDOFile &RAYFile, std::ofstream &ARRFile, DirectOFile &SHDFile)
{
    real atten;
    float freq;
    std::string PlotType;
    
    switch(Beam->RunType[0]){
    case 'R':
    case 'E':
        // Ray trace or Eigenrays
        RAYFile.open(FileRoot + ".ray");
        RAYFile << Title << '\n';
        RAYFile << freqinfo->freq0 << '\n';
        RAYFile << Pos->NSx << Pos->NSy << Pos->NSz << '\n';
        RAYFile << Angles->Nalpha << Angles->Nbeta << '\n';
        RAYFile << Bdry->Top.hs.Depth << '\n';
        RAYFile << Bdry->Bot.hs.Depth << '\n';
        
        RAYFile << (ThreeD ? "xyz" : "rz") << '\n';
        break;
    case 'A':
        // arrival file in ascii format
        ARRFile.open(FileRoot + ".arr");
        
        ARRFile << (ThreeD ? "'3D'" : "'2D'") << "\n";
        
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
        ARRFile.open(FileRoot + ".arr");
        
        ARRFile.write((ThreeD ? "'3D'" : "'2D'"), 4);
        
        freq = freqinfo->freq0; ARRFile.write((char*)&freq, 4);
        
        // write source locations
        
        #define RS sizeof(real)
        if(ThreeD){
            ARRFile.write((char*)&Pos->NSx, 4); for(int32_t i=0; i<Pos->NSx; ++i) ARRFile.write((char*)&Pos->Sx[i], RS);
            ARRFile.write((char*)&Pos->NSy, 4); for(int32_t i=0; i<Pos->NSy; ++i) ARRFile.write((char*)&Pos->Sy[i], RS);
        }
        ARRFile.write((char*)&Pos->NSz, 4); for(int32_t i=0; i<Pos->NSz; ++i) ARRFile.write((char*)&Pos->Sz[i], RS);
        
        // write receiver locations
        ARRFile.write((char*)&Pos->NRz, 4); for(int32_t i=0; i<Pos->NRz; ++i) ARRFile.write((char*)&Pos->Rz[i], RS);
        ARRFile.write((char*)&Pos->NRr, 4); for(int32_t i=0; i<Pos->NRr; ++i) ARRFile.write((char*)&Pos->Rr[i], RS);
        if(ThreeD){
            ARRFile.write((char*)&Pos->Ntheta, 4); for(int32_t i=0; i<Pos->Ntheta; ++i) ARRFile.write((char*)&Pos->theta[i], RS);
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
        WriteHeader(SHDFile, FileRoot + ".shd", Title, atten, PlotType, Pos, freqinfo);
    }
}
