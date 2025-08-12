#include "save_load.h"
#include "fs_compat.h"
#include <fstream>
#include <bits/stdc++.h>
using namespace std;
bool saveProject(const std::string& folderPath, const SCircuit& ckt) {
    cout << "1" << endl;
//    if(!make_dirs(folderPath)) return false;
    cout << "2" << endl;
    std::ofstream os(folderPath);
    cout << "3" << endl;
    if(!os) return false;
    cout << "4" << endl;
    cereal::JSONOutputArchive ar(os);
    ar(CEREAL_NVP(ckt));
    return true;
}

bool loadProject(const std::string& filePath, SCircuit& out){
    std::ifstream is(filePath);
    if(!is) return false;
    cereal::JSONInputArchive ar(is);
    ar(cereal::make_nvp("ckt", out));   // چون حین save کلید "ckt" بود
    return true;
}

