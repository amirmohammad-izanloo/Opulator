#pragma once
#include <string>

#if defined(_WIN32)
#include <direct.h>   // _mkdir
#include <sys/stat.h> // _stat
inline bool dir_exists(const std::string& p){
    struct _stat st{};
    return _stat(p.c_str(), &st)==0 && (st.st_mode & _S_IFDIR);
}
inline int make_dir(const char* p){ return _mkdir(p); }
#else
#include <sys/stat.h>
  #include <sys/types.h>
  #include <unistd.h>
  inline bool dir_exists(const std::string& p){
    struct stat st{};
    return stat(p.c_str(), &st)==0 && S_ISDIR(st.st_mode);
  }
  inline int make_dir(const char* p){ return mkdir(p, 0755); }
#endif

// ساخت پوشه‌ها به‌صورت بازگشتی (C++14 خالص)
inline bool make_dirs(const std::string& path){
    if(path.empty()) return false;
    std::string acc;
    for(char c : path){
        acc.push_back(c);
        if(c=='/' || c=='\\'){
            if(!dir_exists(acc)) {
                if(make_dir(acc.c_str())!=0 && !dir_exists(acc)) return false;
            }
        }
    }
    if(!dir_exists(acc)) {
        if(make_dir(acc.c_str())!=0 && !dir_exists(acc)) return false;
    }
    return true;
}
