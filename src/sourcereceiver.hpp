#pragma once
#include "common.hpp"
#include "subtab.hpp"

struct Position {
    int32_t NSx = 1, NSy = 1, NSz, NRz, NRr, Ntheta; // number of x, y, z, r, theta coordinates
    real Delta_r, Delta_theta;
    int32_t *iSz, *iRz;
    real *Sx, *Sy, *Sz; // Source x, y, z coordinates
    real *Rr, *Rz, *ws, *wr; // Receiver r, z coordinates and weights for interpolation
    real *theta; // Receiver bearings
};

/**
 * Optionally reads a vector of source frequencies for a broadband run
 * If the broadband option is not selected, then the input freq (a scalar) is stored in the frequency vector
 * 
 * freq0: Nominal or carrier frequency
 */
inline void ReadfreqVec(real freq0, char BroadbandOption, 
    LDIFile &ENVFile, std::ostream &PRTFile,
    int32_t &Nfreq, real *&freqVec)
{
    Nfreq = 1;
    
    if(BroadbandOption == 'B'){
        ENVFile.List(); ENVFile.Read(Nfreq);
        PRTFile << "__________________________________________________________________________\n\n\n";
        PRTFile << "Number of frequencies = " << Nfreq << "\n";
        if(Nfreq <= 0){
            std::cout << "Number of frequencies must be positive\n";
            std::abort();
        }
    }
    
    if(freqVec != nullptr) deallocate(freqVec);
    freqVec = allocate<real>(std::max(3, Nfreq));
    
    if(BroadbandOption == 'B'){
        PRTFile << "Frequencies (Hz)\n";
        freqVec[2] = RC(-999.9);
        ENVFile.List(); ENVFile.Read(freqVec, Nfreq);
        SubTab(freqVec, Nfreq);
        EchoVector(freqVec, Nfreq, PRTFile);
    }else{
        freqVec[0] = freq0;
    }
}

/**
 * Read a vector x
 * Description is something like 'receiver ranges'
 * Units       is something like 'km'
 */
inline void ReadVector(int32_t &Nx, real *&x, std::string Description, 
    std::string Units, LDIFile &ENVFile, std::ostream &PRTFile)
{
    PRTFile << "\n__________________________________________________________________________\n\n";
    ENVFile.List(); ENVFile.Read(Nx);
    PRTFile << "Number of " << Description << " = " << Nx << "\n";
    
    if(Nx <= 0){
        std::cout << "ReadVector: Number of " << Description << " must be positive\n";
        std::abort();
    }
    
    if(x != nullptr) deallocate(x);
    x = allocate<real>(Nx);
    
    PRTFile << Description << " (" << Units << ")\n";
    x[2] = RC(-999.9);
    ENVFile.List(); ENVFile.Read(x, Nx);
    
    SubTab(x, Nx);
    Sort(x, Nx);
    EchoVector(x, Nx, PRTFile);
    
    PRTFile << "\n";
    
    // Vectors in km should be converted to m for internal use
    Units.trim();
    if(Units.length() >= 2){
        if(Units.substr(2) == "km") for(int32_t i=0; i<Nx; ++i) x[i] *= RC(1000.0);
    }
}

/**
 * Read source x-y coordinates
 * 
 * ThreeD: flag indicating whether this is a 3D run
 */
inline void ReadSxSy(bool ThreeD, LDIFile &ENVFile, std::ostream &PRTFile,
    Position &Pos)
{
    if(ThreeD){
        ReadVector(Pos.NSx, Pos.Sx, "source   x-coordinates, Sx", "km", ENVFile, PRTFile);
        ReadVector(Pos.NSy, Pos.Sy, "source   y-coordinates, Sy", "km", ENVFile, PRTFile);
    }else{
        Pos.Sx = allocate<real>(1); Pos.Sy = allocate<real>(1);
        Pos.Sx[0] = RC(0.0);
        Pos.Sy[0] = RC(0.0);
    }
}

/**
 * Reads source and receiver z-coordinates (depths)
 * zMin, zMax: limits for those depths; 
 *     sources and receivers are shifted to be within those limits
 */
inline void ReadSzRz(real zMin, real zMax, LDIFile &ENVFile, std::ostream &PRTFile,
    Position &Pos)
{
    //bool monotonic; //LP: monotonic is a function, this is a name clash
    
    ReadVector(Pos.NSz, Pos.Sz, "Source   depths, Sz", "m", ENVFile, PRTFile);
    ReadVector(Pos.NRz, Pos.Rz, "Receiver depths, Rz", "m", ENVFile, PRTFile);
    
    if(Pos.ws != nullptr) deallocate(Pos.ws); deallocate(Pos.iSz);
    Pos.ws = allocate<real>(Pos.NSz); Pos.iSz = allocate<int32_t>(Pos.NSz);
    
    if(Pos.wr != nullptr) deallocate(Pos.wr); deallocate(Pos.iRz);
    Pos.Wr = allocate<real>(Pos.NRz); Pos.iRz = allocate<int32_t>(Pos.NRz);
    
    // *** Check for Sz/Rz in upper or lower halfspace ***
    
    bool topbdry = false, botbdry = false;
    for(int32_t i=0; i<Pos.NSz; ++i){
        if(Pos.Sz[i] < zMin){
            topbdry = true;
            Pos.Sz[i] = zMin;
        }
        if(Pos.Sz[i] > zMax){
            botbdry = true;
            Pos.Sz[i] = zMax;
        }
    }
    if(topbdry) PRTFile << "Warning in ReadSzRz : Source above or too near the top bdry has been moved down\n";
    if(botbdry) PRTFile << "Warning in ReadSzRz : Source below or too near the bottom bdry has been moved up\n";
    
    topbdry = false; botbdry = false;
    for(int32_t i=0; i<Pos.NRz; ++i){
        if(Pos.Rz[i] < zMin){
            topbdry = true;
            Pos.Rz[i] = zMin;
        }
        if(Pos.Rz[i] > zMax){
            botbdry = true;
            Pos.Rz[i] = zMax;
        }
    }
    if(topbdry) PRTFile << "Warning in ReadSzRz : Receiver above or too near the top bdry has been moved down\n";
    if(botbdry) PRTFile << "Warning in ReadSzRz : Receiver below or too near the bottom bdry has been moved up\n";
    
    /*
    if(!monotonic(Pos.sz, Pos.NSz)){
        std::cout << "SzRzRMod: Source depths are not monotonically increasing\n";
        std::abort();
    }
    if(!monotonic(Pos.rz, Pos.NRz)){
        std::cout << "SzRzRMod: Receiver depths are not monotonically increasing\n";
        std::abort();
    }
    */
}

inline void ReadRcvrRanges(LDIFile &ENVFile, std::ostream &PRTFile,
    Position &Pos)
{
    ReadVector(Pos.NRr, Pos.Rr, "Receiver ranges, Rr", "km", ENVFile, PRTFile);
    
    // calculate range spacing
    Pos.delta_r = RC(0.0);
    if(Pos.NRr != 1) Pos.delta_r = Pos.Rr[Pos.NRr-1] - Pos.Rr[Pos.NRr-2];
    
    if(!monotonic(Pos.rr, Pos.NRr)){
        std::cout << "ReadRcvrRanges: Receiver ranges are not monotonically increasing\n";
        std::abort();
    }
}

inline void ReadRcvrBearings(LDIFile &ENVFile, std::ostream &PRTFile,
    Position &Pos)
{
    ReadVector(Pos.Ntheta, Pos.theta, "receiver bearings, theta", "degrees");
    
    CheckFix360Sweep(Pos.theta, Pos.Ntheta);
    
    // calculate angular spacing
    Pos.Delta_theta = RC(0.0);
    if(Pos.Ntheta != 1) Pos.Delta_theta = Pos.theta[Pos.Ntheta-1] - Pos.theta[Pos.Ntheta-2];
    
    if(!monotonic(Pos.theta, Pos.Ntheta)){
        std::cout << "ReadRcvrBearings: Receiver bearings are not monotonically increasing\n";
    }
}
