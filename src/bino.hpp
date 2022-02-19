#pragma once

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common.hpp!"
#endif

/**
 * C++ emulation of FORTRAN direct output (binary). Uses a global record length
 * which is not directly encoded in the file (usually user-encoded as the first
 * word).
 */
class DirectOFile {
public:
    DirectOFile() : recl(-1), bytesWrittenThisRecord(777777777) {}
    ~DirectOFile() { if(ostr.is_open()) ostr.close(); }
    
    /**
     * LRecl: record length in bytes
     */
    void open(const std::string &path, int32_t LRecl){
        recl = LRecl;
        ostr.open(path);
    }
    
    bool good(){
        return ostr.good() && ostr.is_open();
    }
    
    void rec(int32_t r){
        ostr.seekp(r * recl);
        bytesWrittenThisRecord = 0;
    }
    
    #define DOFWRITE(d, data, bytes) d.write(__FILE__, __LINE__, data, bytes)
    void write(const char *file, int fline, const void *data, int32_t bytes){
        if(bytesWrittenThisRecord + bytes > recl){
            std::cout << file << ":" << fline << ": DirectOFile overflow, " 
                << bytesWrittenThisRecord << " bytes already written, rec size "
                << recl << ", tried to write " << bytes << " more\n";
            std::abort();
        }
        ostr.write((const char*)data, bytes);
        bytesWrittenThisRecord += bytes;
    }
    void write(const char *file, int fline, const std::string &str, int32_t bytes){
        write(file, fline, str.data(), math::min(bytes, (int32_t)str.size()));
        for(int32_t b=str.size(); b<bytes; ++b){
            ostr.put(' ');
            ++bytesWrittenThisRecord;
        }
    }
    #define DOFWRITEV(d, data) d.write(__FILE__, __LINE__, data)
    template<typename T> void write(const char *file, int fline, T v){
        write(file, fline, &v, sizeof(T));
    }
    
private:
    std::ofstream ostr;
    int32_t recl;
    int32_t bytesWrittenThisRecord;
};

/**
 * C++ emulation of FORTRAN unformatted output (binary). Each FORTRAN write
 * statement defines a record, and the record length is encoded as an int32_t
 * at the beginning and end of each record.
 */
class UnformattedOFile {
public:
    UnformattedOFile() : recstart(-1), recl(-1) {}
    ~UnformattedOFile() {
        if(ostr.is_open()){
            FinishRecord();
            ostr.close();
        }
    }
    
    void open(const std::string &path) { ostr.open(path); }
    
    bool good() { return ostr.good() && ostr.is_open(); }
    
    void rec(){
        FinishRecord();
        int32_t zero = 0;
        ostr.write((const char*)&zero, 4);
    }
    
    template<typename T> void write(T v){
        if(recstart < 0){
            std::cout << "Missing record in UnformattedOFile!\n";
            bail();
        }
        ostr.write((const char*)&v, sizeof(T));
        recl += sizeof(T);
    }
    
    template<typename T> void write(T *arr, size_t n){
        if(recstart < 0){
            std::cout << "Missing record in UnformattedOFile!\n";
            bail();
        }
        for(size_t i=0; i<n; ++i) ostr.write((const char*)&arr[i], sizeof(T));
        recl += n * sizeof(T);
    }
    
private:
    void FinishRecord() {
        if(recstart < 0 || recl < 0){
            recstart = 0;
            recl = 0;
            return;
        }
        ostr.seekp(recstart);
        ostr.write((const char*)&recl, 4);
        ostr.seekp(recstart+4+recl);
        ostr.write((const char*)&recl, 4);
        recstart += recl + 8;
        recl = 0;
    }
    
    std::ofstream ostr;
    int32_t recstart;
    int32_t recl;
};
