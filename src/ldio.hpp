#pragma once
#include "common.hpp"

#include <fstream>
#include <string>
#include <cctype>

class LDIFile {
public:
    LDIFile(const std::string &filename) : _filename(filename), lastitemcount(0), line(0){
        f.open(filename);
        if(!f.good()) Error("Failed to open file");
        ++line;
        //automatically closed when ifstream destroyed
    }
    
    
    
private:
    void Error(std::string msg){
        std::cout << _filename << ":" << line << ": " << msg << "\n";
        std::abort();
    }
    std::string GetNextItem(){
        if(lastitemcount > 0){
            --lastitemcount;
            return lastitem;
        }
        //Whitespace before start of item
        while(isspace(f.peek()){
            int c = f.get();
            if(c == '\n') ++line;
        }
        if(f.peek() == ','){
            f.get();
            return ""; //Null item
        }
        lastitem = "";
        int quotemode = 0;
        while(true){
            int c = f.get();
            if(c == '"'){
                if(!quotemode){
                    quotemode = 1;
                }else if(quotemode == 1){
                    break; //End of string
                }else{
                    lastitem += c; //add opposite quote to string
                }
            }else if(c == '\''){
                if(!quotemode){
                    quotemode = 2;
                }else if(quotemode == 2){
                    break; //End of string
                }else{
                    lastitem += c; //add opposite quote to string
                }
            }else if(isspace(c)){
                if(c == '\n') ++line;
                break; //End of item
            }else if(c == '*'){
                if(lastitemcount != 0) Error("Can't have nested repetitions");
                if(!isInt(lastitem, false)) Error("Invalid repetition count");
                lastitemcount = std::stoul(lastitem);
                if(lastitemcount == 0) Error("Repetition count can't be 0");
                lastitem = "";
            }else{
                lastitem += c;
            }
        }
        //TODO
        if(lastitemcount > 0) --lastitemcount;
        return lastitem;
    }
    
    std::string _filename;
    std::ifstream f;
    std::string lastitem;
    uint32_t lastitemcount;
    uint32_t line;
}; 
