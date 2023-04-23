/*
bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP(3D) underwater acoustics simulator
Copyright (C) 2021-2023 The Regents of the University of California
Marine Physical Lab at Scripps Oceanography, c/o Jules Jaffe, jjaffe@ucsd.edu
Based on BELLHOP / BELLHOP3D, which is Copyright (C) 1983-2022 Michael B. Porter

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
#error "Must be included from common_setup.hpp!"
#endif

namespace bhc {

/**
 * C++ emulation of FORTRAN direct output (binary). Uses a global record length
 * which is not directly encoded in the file (usually user-encoded as the first
 * word).
 */
class DirectOFile {
public:
    DirectOFile(bhcInternal *internal)
        : _internal(internal), recl(0), record(777777777),
          bytesWrittenThisRecord(777777777), highestRecord(0),
          bytesWrittenHighestRecord(0)
    {}
    ~DirectOFile()
    {
        if(!ostr.is_open()) return;
        // End of highest record may not have been filled up, causing the file
        // length to be too short, if only part of it was written
        if(bytesWrittenHighestRecord < recl) {
            ostr.seekp((highestRecord + 1) * recl - 1);
            ostr.put('\0');
        }
        ostr.close();
    }

    /**
     * LRecl: record length in bytes
     */
    void open(const std::string &path, size_t LRecl)
    {
        recl = LRecl;
        ostr.open(path, std::ios::binary);
    }

    bool good() { return ostr.good() && ostr.is_open(); }

    void rec(size_t r)
    {
        if(r >= highestRecord) {
            highestRecord             = r;
            bytesWrittenHighestRecord = 0;
        }
        record = r;
        ostr.seekp(r * recl);
        bytesWrittenThisRecord = 0;
    }

#define DOFWRITE(d, data, bytes) d.write(__FILE__, __LINE__, data, bytes)
    void write(const char *file, int fline, const void *data, size_t bytes)
    {
        checkAndIncrement(file, fline, bytes);
        ostr.write((const char *)data, bytes);
    }
    void write(const char *file, int fline, const std::string &str, size_t bytes)
    {
        checkAndIncrement(file, fline, bytes);
        ostr.write(str.data(), bhc::min(bytes, str.size()));
        for(size_t b = str.size(); b < bytes; ++b) ostr.put(' ');
    }
#define DOFWRITEV(d, data) d.write(__FILE__, __LINE__, data)
    template<typename T> void write(const char *file, int fline, T v)
    {
        write(file, fline, &v, sizeof(T));
    }

private:
    bhcInternal *_internal;
    std::ofstream ostr;
    size_t recl;
    size_t record;
    size_t bytesWrittenThisRecord;
    size_t highestRecord;
    size_t bytesWrittenHighestRecord;

    void checkAndIncrement(const char *file, int fline, size_t bytes)
    {
        if(bytesWrittenThisRecord + bytes > recl) {
            ExternalError(
                _internal,
                "%s:%d: DirectOFile overflow, %" PRIuMAX
                " bytes already written, rec size %" PRIuMAX ", tried to write %" PRIuMAX
                " more",
                file, fline, bytesWrittenThisRecord, recl, bytes);
        }
        bytesWrittenThisRecord += bytes;
        if(record == highestRecord) {
            bytesWrittenHighestRecord = bytesWrittenThisRecord;
        }
    }
};

class DirectIFile {
public:
    DirectIFile(bhcInternal *internal) : _internal(internal), recl(0) {}
    ~DirectIFile()
    {
        if(!istr.is_open()) return;
        istr.close();
    }

    void open(const std::string &path)
    {
        istr.open(path, std::ios::binary);
        if(!(istr.good() && istr.is_open())) {
            ExternalError(_internal, "Failed to open DirectIFile %s", path.c_str());
        }
        int32_t temp;
        istr.read((char *)&temp, 4);
        recl = temp * 4;
        istr.seekg(0, istr.end);
        fileLen = istr.tellg();
        if((fileLen % recl) != 0) {
            ExternalError(
                _internal,
                "DirectIFile %s file length %" PRIuMAX
                " is not multiple of record length %" PRIuMAX,
                path.c_str(), fileLen, recl);
        }
    }

#define DIFREC(d, r) d.rec(__FILE__, __LINE__, r)
    void rec(const char *file, int fline, size_t r)
    {
        size_t a = r * recl;
        if(a >= fileLen) {
            ExternalError(
                _internal,
                "%s:%d: DirectIFile record %" PRIuMAX
                " out of bounds, record length is %" PRIuMAX ", file length is %" PRIuMAX,
                file, fline, r, recl, fileLen);
        }
        record = r;
        istr.seekg(a);
        bytesReadThisRecord = 0;
    }

#define DIFREAD(d, data, bytes) d.read(__FILE__, __LINE__, data, bytes)
    void read(const char *file, int fline, void *data, size_t bytes)
    {
        checkAndIncrement(file, fline, bytes);
        istr.read((char *)data, bytes);
    }

#define DIFREADS(d, bytes) d.readstring(__FILE__, __LINE__, bytes)
    std::string readstring(const char *file, int fline, size_t bytes)
    {
        char *s = new char[bytes];
        read(file, fline, s, bytes);
        std::string ret(s, bytes);
        delete[] s;
        return ret;
    }

#define DIFREADV(d, data) d.read(__FILE__, __LINE__, data)
    template<typename T> void read(const char *file, int fline, T &v)
    {
        read(file, fline, &v, sizeof(T));
    }

#define DIFSKIP(d, bytes) d.skip(__FILE__, __LINE__, bytes)
    void skip(const char *file, int fline, size_t bytes)
    {
        checkAndIncrement(file, fline, bytes);
        istr.seekg((int)bytes, istr.cur);
    }

private:
    bhcInternal *_internal;
    std::ifstream istr;
    size_t recl;
    size_t record;
    size_t bytesReadThisRecord;
    size_t fileLen;

    void checkAndIncrement(const char *file, int fline, size_t bytes)
    {
        if(bytesReadThisRecord + bytes > recl) {
            ExternalError(
                _internal,
                "%s:%d: DirectIFile overflow, %" PRIuMAX
                " bytes already read, rec size %" PRIuMAX ", tried to read %" PRIuMAX
                " more",
                file, fline, bytesReadThisRecord, recl, bytes);
        }
        bytesReadThisRecord += bytes;
    }
};

} // namespace bhc
