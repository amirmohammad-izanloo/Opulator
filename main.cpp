#include <iostream>
#include <string>
#include <vector>
#include <regex>
using namespace std;
//--------------------- parse commands ---------------------
vector<string> parseCommands(const string& command) {
    vector<string> params;
    vector<regex> patterns = {
        regex(R"(add\s+R([^<]+)\s+([^<]+)\s+([^<]+)\s+([^<]+))"),
        regex(R"(delete\s+R([^<]+))"),
        regex(R"(add\s+C([^<]+)\s+([^<]+)\s+([^<]+)\s+([^<]+))"),
        regex(R"(delete\s+C([^<]+))"),
        regex(R"(add\s+L([^<]+)\s+([^<]+)\s+([^<]+)\s+([^<]+))"),
        regex(R"(delete\s+L([^<]+))"),
        regex(R"(add\s+D([^<]+)\s+([^<]+)\s+([^<]+)\s+([^<]+))"),
        regex(R"(delete\s+D([^<]+))"),
        regex(R"(add\s+GND\s+([^<]+))"),
        regex(R"(delete\s+GND\s+([^<]+))"),
        regex(R"(list\s+([^<]+))"),
        regex(R"(rename\s+node\s+([^<]+)\s+([^<]+))")
    };
    smatch match;
    for (const auto& pat : patterns) {
        if (regex_match(command, match, pat)) {
            for (size_t i = 1; i < match.size(); ++i)
                params.push_back(match[i]);
            return params;
        }
    }
    return {};
}
//-------------------- parse number ---------------------------
double parseNumber(string input) {
    double result = 0;
    char last = input[input.length() - 1];
    string numberPart;
    if (last == 'k' || last == 'K') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) * 1000;
    } else if (last == 'm' || last == 'M') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) / 1000;
    } else {
        result = stod(input);
    }
    return result;
}
//------------------- node ----------------------
struct Node {
    string name;
    int voltage;
    Node(): name(""), voltage(0) {}
};
//------------------- classes -------------------
class Element {
protected:
    string name;
    Node n1, n2;
    double value;
    string unit;
public:
    Element(const string &elemName, const Node &node1, const Node &node2, double val, const string &unitType)
        : name(elemName), n1(node1), n2(node2), value(val), unit(unitType) {}

    virtual ~Element() {}
};

class Resistor : public Element {
    double conductance;
public:
    Resistor(const string &elemName, const Node &node1, const Node &node2, double resistance)
        : Element(elemName, node1, node2, resistance, "Ohm")
    {
        conductance = (resistance != 0) ? (1.0 / resistance) : 0;
    }

};

class Capacitor : public Element {
public:
    Capacitor(const string &elemName, const Node &node1, const Node &node2, double capacitance)
        : Element(elemName, node1, node2, capacitance, "Farad") {}

};

class Inductor : public Element {
public:
    Inductor(const string &elemName, const Node &node1, const Node &node2, double inductance)
        : Element(elemName, node1, node2, inductance, "Henry") {}

};

class VoltageSource : public Element {
public:
    VoltageSource(const string &elemName, const Node &node1, const Node &node2, double voltage)
        : Element(elemName, node1, node2, voltage, "Volt") {
    };
};

class CurrentSource : public Element {
public:
    CurrentSource(const string &elemName, const Node &node1, const Node &node2, double current)
        : Element(elemName, node1, node2, current, "Ampere") {}
};
//-----------------------------------------------
int main() {

    cout << "Hello Izan" << endl;
    return 0;
}
