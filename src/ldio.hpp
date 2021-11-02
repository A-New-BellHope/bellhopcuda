#pragma once
#include "common.hpp"

/**
 * C++ emulation of FORTRAN list-directed input.
 * To use:
 * LDIFile YourFile("filename");
 * //List() starts a new list (single READ line). This is needed because the
 * //slash in an input file terminates assignment for the rest of the list,
 * //so we need to be able to distinguish list boundaries.
 * YourFile.List(); YourFile.Read(somestring); YourFile.read(somecpx); //etc.
 */
class LDIFile {
public:
    LDIFile(const std::string &filename, bool abort_on_error = true) 
        : _filename(filename), _abort_on_error(abort_on_error),
        lastitemcount(0), line(0), isafterslash(false)
    {
        f.open(filename);
        if(!f.good()) Error("Failed to open file");
        ++line;
        //automatically closed when ifstream destroyed
    }
    
    bool Good() { return f.good(); }
    
    void List(){
        isafterslash = false;
    }
    
    #define LDIFILE_READPREFIX() \
        if(f.eof() && !isafterslash) Error("End of file"); \
        std::string s = GetNextItem(); \
        if(s == nullitem) return; \
        do{} while(false)
        
    void Read(std::string &v){
        LDIFILE_READPREFIX();
        v = s;
    }
    void Read(int32_t &v){
        LDIFILE_READPREFIX();
        if(!isInt(s, true)) Error("String " + s + " is not an integer");
        v = std::stoi(s);
    }
    void Read(uint32_t &v){
        LDIFILE_READPREFIX();
        if(!isInt(s, false)) Error("String " + s + " is not an unsigned integer");
        v = std::stoul(s);
    }
    void Read(real &v){
        LDIFILE_READPREFIX();
        if(!isReal(s)) Error("String " + s + " is not a real number");
        v = strtod(s.c_str(), nullptr);
    }
    void Read(vec2 &v){
        LDIFILE_READPREFIX();
        if(!isReal(s)) Error("String " + s + " is not a real number");
        v.x = strtod(s.c_str(), nullptr);
        if(f.eof() && !isafterslash) Error("End of file");
        s = GetNextItem();
        if(s == nullitem) Error("Only specified part of a vec2!");
        if(!isReal(s)) Error("String " + s + " is not a real number");
        v.y = strtod(s.c_str(), nullptr);
    }
    void Read(cpx &v){
        LDIFILE_READPREFIX();
        if(s[0] != '(' || s.back() != ')') 
            Error("String " + s + " is not a complex number (not in parentheses)");
        size_t commapos = s.find(',');
        if(commapos == std::string::npos || s.find(',', commapos+1) != std::string::npos)
            Error("String " + s + " is not a complex number (not exactly one comma)");
        std::string sr, si;
        sr = s.substr(1, commapos-1);
        si = s.substr(commapos+1, s.length()-commapos-2);
        if(!isReal(sr) || !isReal(si))
            Error("String " + s + " is not a complex number (components not real numbers)");
        real vr = strtod(sr.c_str(), nullptr);
        real vi = strtod(si.c_str(), nullptr);
        v = cpx(vr, vi);
    }
    void Read(char *v, size_t nc){
        LDIFILE_READPREFIX();
        if(s.length() > nc) Error("Max string length " + std::to_string(nc) 
            + ", got " + s);
        std::strncpy(v, s.c_str(), s.length());
        if(s.length() < nc) memset(v+s.length(), ' ', nc-s.length());
    }
    void Read(real *arr, size_t count){
        for(size_t i=0; i<count; ++i){
            Read(arr[i]);
        }
    }
    
private:
    void Error(std::string msg){
        std::cout << _filename << ":" << line << ": " << msg << "\n";
        std::cout << "Last token is: \"" << lastitem << "\"\n";
        if(_abort_on_error) std::abort();
    }
    void IgnoreRestOfLine(){
        while(!f.eof() && f.peek() != '\n') f.get();
        if(!f.eof()) f.get(); //get the \n
        ++line;
    }
    void GotSlash(){
        isafterslash = true;
        IgnoreRestOfLine();
    }
    std::string GetNextItem(){
        if(lastitemcount > 0){
            --lastitemcount;
            //std::cout << "-- lastitemcount, returning " << lastitem << "\n";
            return lastitem;
        }
        if(isafterslash){
            //std::cout << "-- isafterslash, returning null\n";
            return nullitem;
        }
        //Whitespace before start of item
        while(!f.eof()){
            int c = f.peek();
            if(c == '!'){
                IgnoreRestOfLine();
                continue;
            }else if(isspace(c)){
                if(c == '\n') ++line;
            }else{
                break;
            }
            f.get();
        }
        if(f.eof()) return nullitem;
        if(f.peek() == ','){
            f.get();
            //std::cout << "-- empty comma, returning null\n";
            return nullitem;
        }
        //Main item
        lastitem = "";
        int quotemode = 0;
        while(!f.eof()){
            int c = f.get();
            if(quotemode == 1){
                if(c == '"'){
                    quotemode = -1;
                    break;
                }else if(c == '\n'){
                    ++line;
                    break;
                }else{
                    lastitem += c;
                }
            }else if(quotemode == 2){
                if(c == '\''){
                    quotemode = -1;
                    break;
                }else if(c == '\n'){
                    ++line;
                    break;
                }else{
                    lastitem += c;
                }
            }else{
                if(c == '"'){
                    quotemode = 1;
                }else if(c == '\''){
                    quotemode = 2;
                }else if(c == '('){
                    if(!quotemode){
                        quotemode = 3;
                    }
                    lastitem += c;
                }else if(c == ')'){
                    if(quotemode == 3){
                        quotemode = -1;
                        lastitem += c;
                        break;
                    }
                    lastitem += c;
                }else if(isspace(c)){
                    if(c == '\n') ++line;
                    break;
                }else if(c == ',' && quotemode != 3){
                    break;
                }else if(c == '*'){
                    if(lastitemcount != 0) Error("Can't have nested repetitions");
                    if(!isInt(lastitem, false)) Error("Invalid repetition count");
                    lastitemcount = std::stoul(lastitem);
                    if(lastitemcount == 0) Error("Repetition count can't be 0");
                    lastitem = "";
                }else if(c == '/'){
                    GotSlash();
                    break;
                }else if(c == '!'){
                    IgnoreRestOfLine();
                    break;
                }else{
                    lastitem += c;
                }
            }
        }
        if(quotemode > 0) Error("Quotes or parentheses not closed");
        if(f.eof()){
            //std::cout << "-- eof, returning " << lastitem << "\n";
            return lastitem;
        }
        if(quotemode < 0){
            int c = f.peek();
            if(!isspace(c) && c != ',') Error(std::string("Invalid character '")
                + (char)c + std::string("' after end of quoted string"));
        }
        if(isafterslash){
            //std::cout << "-- new isafterslash, returning " << lastitem << "\n";
            return lastitem;
        }
        //Whitespace and comma after item
        bool hadcomma = false;
        while(!f.eof()){
            int c = f.peek();
            if(c == '!'){
                IgnoreRestOfLine();
                continue;
            }else if(isspace(c)){
                if(c == '\n') ++line;
            }else if(c == ','){
                if(!hadcomma){
                    hadcomma = true;
                }else{
                    break;
                }
            }else if(c == '/'){
                GotSlash();
                break;
            }else{
                break;
            }
            f.get();
        }
        //Finally
        if(lastitemcount > 0) --lastitemcount;
        //std::cout << "-- normal returning " << lastitem << "\n";
        return lastitem;
    }
    
    std::string _filename;
    bool _abort_on_error;
    std::ifstream f;
    std::string lastitem;
    uint32_t lastitemcount;
    uint32_t line;
    bool isafterslash;
    
    //This string is not possible to represent, so we use it to indicate null
    //(empty string is separate and valid)
    static constexpr const char *const nullitem = "\"'";
}; 
