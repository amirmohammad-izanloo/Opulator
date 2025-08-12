#pragma once
#include <vector>
#include <string>

// cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

struct SNode {
    int id{}; std::string name; int x{0}, y{0};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(id), CEREAL_NVP(name), CEREAL_NVP(x), CEREAL_NVP(y));
    }
};

struct Point {
    int x{0}, y{0};
    template<class Ar> void serialize(Ar& ar) {
        ar(CEREAL_NVP(x), CEREAL_NVP(y));
    }
};

struct SWire {
    int id{}; std::vector<Point> points; // polyline
    template<class Ar> void serialize(Ar& ar){ ar(CEREAL_NVP(id), CEREAL_NVP(points)); }
};

struct SElement {
    std::string kind, name;
    int n1{-1}, n2{-1};
    double value{0.0};
    int x{0}, y{0};
    double rotation{0.0}; bool mirror{false};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(kind), CEREAL_NVP(name),
           CEREAL_NVP(n1), CEREAL_NVP(n2), CEREAL_NVP(value),
           CEREAL_NVP(x), CEREAL_NVP(y), CEREAL_NVP(rotation), CEREAL_NVP(mirror));
    }
};

struct TransientCfg {
    double t_stop{1e-3}, t0{0}, t_step{1e-6};
    template<class Ar> void serialize(Ar& ar){ ar(CEREAL_NVP(t_stop), CEREAL_NVP(t0), CEREAL_NVP(t_step)); }
};

struct ACSweepCfg {
    std::string type{"Linear"}; double w_start{10}, w_stop{1e6}; int points{200};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(type), CEREAL_NVP(w_start), CEREAL_NVP(w_stop), CEREAL_NVP(points));
    }
};

struct PhaseSweepCfg {
    double w_base{1000}, phi_start{0}, phi_stop{3.14159}; int points{181};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(w_base), CEREAL_NVP(phi_start), CEREAL_NVP(phi_stop), CEREAL_NVP(points));
    }
};

struct ProjectSettings {
    TransientCfg transient; ACSweepCfg ac; PhaseSweepCfg phase;
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(transient), CEREAL_NVP(ac), CEREAL_NVP(phase));
    }
};

struct SCircuit {
    int file_version{1};
    std::vector<SNode>    Snodes;
    std::vector<SWire>    Swires;
    std::vector<SElement> Selements;
    ProjectSettings       settings;

    template<class Ar> void serialize(Ar& ar){
        ar( cereal::make_nvp("version",  file_version),
            cereal::make_nvp("nodes",    Snodes),
            cereal::make_nvp("wires",    Swires),
            cereal::make_nvp("elements", Selements),
            cereal::make_nvp("settings", settings) );
    }
};

// API
bool saveProject(const std::string& folderPath, const SCircuit& ckt);
bool loadProject(const std::string& folderPath, SCircuit& out);
