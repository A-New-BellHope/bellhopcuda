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

class PrintFileEmu {
public:
    PrintFileEmu(const char *FileRoot, void (*outputCallback)(const char *message))
        : callback(nullptr)
    {
        if(outputCallback == nullptr){
            std::string s(FileRoot);
            ofs.open(s + ".prt");
            if(!ofs.good()){
                std::cout << "Could not open print file: " << FileRoot << ".prt\n";
                std::abort();
            }
            ofs << std::unitbuf;
        }else{
            callback = outputCallback;
        }
    }
    ~PrintFileEmu()
    {
        if(ofs.is_open()) ofs.close();
    }
    
    template<typename T> PrintFileEmu &operator<<(const T &x)
    {
        if(callback != nullptr){
            linebuf << x;
            size_t newlinepos = linebuf.str().find("\n");
            if(newlinepos != std::string::npos){
                std::string msg = linebuf.str().substr(0, newlinepos);
                callback(msg.c_str());
                linebuf.ignore(newlinepos+1);
            }
        }else if(ofs.good()){
            ofs << x;
        }
        return *this;
    }
    
private:
    std::ofstream ofs;
    std::stringstream linebuf;
    void (*callback)(const char *message);
};

}
