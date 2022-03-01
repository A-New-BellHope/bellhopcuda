
/**
 * Write header to disk file
 * LP: of a SHDFile ("binary `shade' file (SHDFIL) that contains calculated
 * pressure fields")
 *
 * FileName: Name of the file (could be a shade file or a Green's function file)
 * Title: Arbitrary title
 * freq0: Nominal frequency [LP: now in freqinfo]
 * Atten: stabilizing attenuation (for wavenumber integration only)
 * PlotType: If "TL", writes only first and last Sx and Sy [LP: never set to
 * "TL" in BELLHOP]
 */
void WriteHeader(DirectOFile &SHDFile, const std::string &FileName, 
    const std::string &Title, float Atten, const std::string &PlotType, 
    const Position *Pos, const FreqInfo *freqinfo)
{
    bool isTL = (PlotType[0] == 'T' && PlotType[1] == 'L');
    
    int32_t LRecl = 84; //4 for LRecl, 80 for Title
    LRecl = math::max(LRecl, 2 * freqinfo->Nfreq * (int32_t)sizeof(freqinfo->freqVec[0]));
    LRecl = math::max(LRecl, Pos->Ntheta * (int32_t)sizeof(Pos->theta[0]));
    if(!isTL){
        LRecl = math::max(LRecl, Pos->NSx * (int32_t)sizeof(Pos->Sx[0]));
        LRecl = math::max(LRecl, Pos->NSy * (int32_t)sizeof(Pos->Sy[0]));
    }
    LRecl = math::max(LRecl, Pos->NSz * (int32_t)sizeof(Pos->Sz[0]));
    LRecl = math::max(LRecl, Pos->NRz * (int32_t)sizeof(Pos->Rz[0]));
    LRecl = math::max(LRecl, Pos->NRr * (int32_t)sizeof(cpxf));
    
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
    LDOFile &RAYFile, DirectOFile &SHDFile)
{
    real atten;
    std::string PlotType;
    
    
        
        
        break;
    case 'A':
        // arrival file in ascii format
        // LP: Moved to WriteArrivals
        break;
    case 'a':
        // arrival file in binary format
        // LP: Moved to WriteArrivals
        break;
    default:
        atten = FL(0.0);
        
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
