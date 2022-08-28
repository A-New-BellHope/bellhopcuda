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
#pragma once

#ifndef _BHC_INCLUDING_COMPONENTS_
#error "Must be included from common.hpp!"
#endif

namespace bhc {

/**
 * C++ emulation of FORTRAN direct output (binary). Uses a global record length
 * which is not directly encoded in the file (usually user-encoded as the first
 * word).
 */
class DirectOFile {
public:
    DirectOFile() : recl(0), bytesWrittenThisRecord(777777777) {}
    ~DirectOFile() { if(ostr.is_open()) ostr.close(); }
    
    /**
     * LRecl: record length in bytes
     */
    void open(const std::string &path, size_t LRecl){
        recl = LRecl;
        ostr.open(path);
    }
    
    bool good(){
        return ostr.good() && ostr.is_open();
    }
    
    void rec(size_t r){
        ostr.seekp(r * recl);
        bytesWrittenThisRecord = 0;
    }
    
    #define DOFWRITE(d, data, bytes) d.write(__FILE__, __LINE__, data, bytes)
    void write(const char *file, int fline, const void *data, size_t bytes){
        if(bytesWrittenThisRecord + bytes > recl){
            GlobalLog("%s:%d: DirectOFile overflow, %lli bytes already written, rec size %lli, tried to write %lli more\n",
                file, fline, bytesWrittenThisRecord, recl, bytes);
            std::abort();
        }
        ostr.write((const char*)data, bytes);
        bytesWrittenThisRecord += bytes;
    }
    void write(const char *file, int fline, const std::string &str, size_t bytes){
        write(file, fline, str.data(), bhc::min(bytes, str.size()));
        for(size_t b=str.size(); b<bytes; ++b){
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
    size_t recl;
    size_t bytesWrittenThisRecord;
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
            GlobalLog("Missing record in UnformattedOFile!\n");
            bail();
        }
        ostr.write((const char*)&v, sizeof(T));
        recl += (int32_t)sizeof(T);
    }
    
    template<typename T> void write(T *arr, size_t n){
        if(recstart < 0){
            GlobalLog("Missing record in UnformattedOFile!\n");
            bail();
        }
        for(size_t i=0; i<n; ++i) ostr.write((const char*)&arr[i], sizeof(T));
        recl += (int32_t)(n * sizeof(T));
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

}
