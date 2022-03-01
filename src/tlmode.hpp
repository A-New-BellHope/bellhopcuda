
void OpenOutputFiles(std::string FileRoot, bool ThreeD, std::string Title,
    const BdryType *Bdry, const Position *Pos, const AnglesStructure *Angles, 
    const FreqInfo *freqinfo, const BeamStructure *Beam,
    LDOFile &RAYFile, DirectOFile &SHDFile);
    
    

    OpenOutputFiles(FileRoot, false, Title, Bdry, Pos, Angles, freqinfo, Beam,
        RAYFile, SHDFile);
