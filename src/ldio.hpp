#pragma once
#include "common.hpp"

/**
 * C++ emulation of FORTRAN list-directed input.
 * To use:
 * LDIFile YourFile("filename");
 * LIST(YourFile); YourFile.Read(somestring); YourFile.read(somecpx); //etc.
 * 
 * List() starts a new list (single READ line). This is needed because ending a
 * list advances to the next line of the input. The converse is not true; a
 * newline does not terminate assignment of values to a list, it will continue
 * reading onto future lines, unless there is a '/'.
 * 
 * LIST_WARNLINE() is for cases when the input variables should all be on the
 * same line of the input file. If this option is used and reading goes onto
 * a new line of the input, a warning is printed.
 */
class LDIFile {
public:
    LDIFile(const std::string &filename, bool abort_on_error = true) 
        : _filename(filename), _abort_on_error(abort_on_error),
        lastitemcount(0), line(0), isafterslash(false), isafternewline(true)
    {
        f.open(filename);
        if(!f.good()) Error("Failed to open file");
        ++line;
        //automatically closed when ifstream destroyed
    }
    
    bool Good() { return f.good(); }
    
#define LIST(ldif) ldif.List(__FILE__, __LINE__)
#define LIST_WARNLINE(ldif) ldif.List(__FILE__, __LINE__, true)
    void List(const char *file, int fline, bool warnline = false){
        codefile = SOURCE_FILENAME(file);
        codeline = fline;
        isafterslash = false;
        if(!isafternewline){
            IgnoreRestOfLine();
        }
        if(warnline){
            _warnline = line;
        }else{
            _warnline = -1;
        }
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
    void Read(char &v){
        LDIFILE_READPREFIX();
        if(s.length() != 1) Error("String " + s + " is not one character");
        v = s[0];
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
    void Read(float &v){
        LDIFILE_READPREFIX();
        if(!isReal(s)) Error("String " + s + " is not a real number");
        v = strtof(s.c_str(), nullptr);
    }
    void Read(double &v){
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
    template<typename REAL> void Read(REAL *arr, size_t count){
        for(size_t i=0; i<count; ++i){
            Read(arr[i]);
        }
    }
    
private:
    void PrintLoc(){
        std::cout << codefile << ":" << codeline << " reading " 
            << _filename << ":" << line  << ": ";
    }
    void Error(std::string msg){
        PrintLoc();
        std::cout << msg << "\n";
        std::cout << "Last token is: \"" << lastitem << "\"\n";
        if(_abort_on_error) std::abort();
    }
    void IgnoreRestOfLine(){
        if(_debug) std::cout << "-- ignoring rest of line\n";
        while(!f.eof() && f.peek() != '\n') f.get();
        if(!f.eof()) f.get(); //get the \n
        ++line;
        isafternewline = true;
    }
    std::string GetNextItem(){
        if(lastitemcount > 0){
            --lastitemcount;
            if(_debug) std::cout << "-- lastitemcount, returning " << lastitem << "\n";
            return lastitem;
        }
        if(isafterslash){
            if(_debug) std::cout << "-- isafterslash, returning null\n";
            return nullitem;
        }
        //Whitespace before start of item
        while(!f.eof()){
            int c = f.peek();
            if(!isspace(c)) break;
            f.get();
            if(c == '\n'){
                ++line;
                isafternewline = true;
            }else{
                isafternewline = false;
            }
        }
        if(f.eof()) return nullitem;
        if(f.peek() == ','){
            f.get();
            isafternewline = false;
            if(_debug) std::cout << "-- empty comma, returning null\n";
            return nullitem;
        }
        //Main item
        if(_warnline >= 0 && _warnline != line){
            PrintLoc();
            std::cout << "Warning: input continues onto next line, likely mistake\n";
            _warnline = line;
        }
        lastitem = "";
        int quotemode = 0;
        while(!f.eof()){
            int c = f.get();
            isafternewline = false;
            if(quotemode == 1){
                if(c == '"'){
                    quotemode = -1;
                    break;
                }else if(c == '\n'){
                    ++line;
                    isafternewline = true;
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
                    isafternewline = true;
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
                    if(c == '\n'){
                        ++line;
                        isafternewline = true;
                    }
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
                    isafterslash = true;
                    break;
                }else{
                    lastitem += c;
                }
            }
        }
        if(quotemode > 0) Error("Quotes or parentheses not closed");
        if(f.eof()){
            if(_debug) std::cout << "-- eof, returning " << lastitem << "\n";
            return lastitem;
        }
        if(quotemode < 0){
            int c = f.peek();
            if(!isspace(c) && c != ',') Error(std::string("Invalid character '")
                + (char)c + std::string("' after end of quoted string"));
        }
        if(isafternewline){
            if(_debug) std::cout << "-- isafternewline, returning " << lastitem << "\n";
            return lastitem;
        }
        if(isafterslash){
            if(_debug) std::cout << "-- new isafterslash, returning " << lastitem << "\n";
            return lastitem;
        }
        //Whitespace and comma after item
        bool hadcomma = false;
        while(!f.eof()){
            int c = f.peek();
            if(isspace(c)){
                f.get();
                if(c != '\n'){
                    isafternewline = false;
                    continue;
                }
                ++line;
                isafternewline = true;
            }else if(c == ','){
                if(!hadcomma){
                    f.get();
                    hadcomma = true;
                    isafternewline = false;
                    continue;
                }
            }else if(c == '/'){
                f.get();
                isafterslash = true;
            }
            break;
        }
        //Finally
        if(lastitemcount > 0) --lastitemcount;
        if(_debug) std::cout << "-- normal returning " << lastitem << "\n";
        return lastitem;
    }
    
    constexpr static bool _debug = false;
    std::string _filename;
    std::string codefile;
    int codeline;
    int _warnline;
    bool _abort_on_error;
    std::ifstream f;
    std::string lastitem;
    uint32_t lastitemcount;
    uint32_t line;
    bool isafterslash, isafternewline;
    
    //This string is not possible to represent, so we use it to indicate null
    //(empty string is separate and valid)
    static constexpr const char *const nullitem = "\"'";
}; 

class LDOFile {
public:
    LDOFile() : iwidth(12), rwidth(24) {}
    ~LDOFile(){
        if(ostr.is_open()) ostr.close();
    }
    
    void open(const std::string &path){
        ostr.open(path);
        ostr << std::setfill(' ');
    }
    bool good(){
        return ostr.good() && ostr.is_open();
    }
    
    LDOFile &operator<<(const char &c){
        ostr << c;
        return *this;
    }
    LDOFile &operator<<(const std::string &s){
        ostr << "'" << s << "'";
        return *this;
    }
    
    void intwidth(int32_t iw) { iwidth = iw; }
    LDOFile &operator<<(const int32_t &i){
        if(iwidth > 0){
            ostr << std::setiosflags(std::ios_base::right)
                << std::setw(iwidth);
        }
        ostr << i;
        return *this;
    }
    
    void realwidth(int32_t rw) { rwidth = rw; }
    LDOFile &operator<<(real r){
        if(rwidth <= 0){
            ostr << r;
            return *this;
        }
        ostr << "  ";
        if(!isfinite(r)){
            ostr << std::setw(rwidth)
                << std::left
                << r;
            return *this;
        }
        bool sci = r != RC(0.0) && (std::abs(r) < RC(0.1) || std::abs(r) >= RC(1.0e6));
        int32_t w = rwidth;
        if(r < RC(0.0)){
            ostr << "-";
            r = -r;
        }else if(sci || r >= RC(1.0) || r == RC(0.0)){
            ostr << " ";
        }
        --w;
        if(sci) --w;
        ostr << std::setprecision(w - 6); //5 for exp, 1 for decimal point
        if(sci){
            ostr << std::setiosflags(std::ios_base::uppercase
                    | std::ios_base::scientific)
                << std::setw(w + 1)
                << std::left
                << r;
        }else{
            ostr << std::setiosflags(std::ios_base::showpoint)
                << std::defaultfloat
                << std::setw(w - 5)
                << r;
            ostr << "     ";
        }
        return *this;
    }
    
    LDOFile &operator<<(const cpx &c){
        ostr << "(";
        this->operator<<(c.real());
        ostr << ",";
        this->operator<<(c.imag());
        ostr << ")";
        return *this;
    }
    LDOFile &operator<<(const vec2 &v){
        this->operator<<(v.x);
        this->operator<<(v.y);
        return *this;
    }
    
    template<typename T> void write(const T *v, int32_t n){
        for(int32_t i=0; i<n; ++i){
            this->operator<<(v[i]);
        }
    }
    template<typename T> void write(const T *v, int32_t m, int32_t n){
        for(int32_t i=0; i<m; ++i){
            for(int32_t j=0; j<n; ++j){
                this->operator<<(v[i*n+j]);
            }
            this->operator<<('\n');
        }
    }
    
private:
    std::ofstream ostr;
    int32_t iwidth, rwidth;
};

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
