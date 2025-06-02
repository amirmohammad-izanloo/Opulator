#include <iostream>
#include <string>
#include <vector>
#include <regex>
using namespace std;


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
    } else if (last == 'u' || last == 'U') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) / 1000000;
    } else {
        result = stod(input);
    }
    return result;
}
//------------------- node ----------------------
class Node {
    string name;
    int voltage;

public:
    Node(const string& nodeName = "", int volt = 0)
            : name(nodeName), voltage(volt) {}

    string getName() const { return name; }
    void setName(const string& newName) { name = newName; }

    int getVoltage() const { return voltage; }
    void setVoltage(int v) { voltage = v; }
};
//------------------- classes -------------------
class Element {
protected:
    string name;
    Node* n1;
    Node* n2;
    double value;
    string unit;
    string type;
public:
    Element(const string &elemName, Node* node1,  Node* node2, double val, const string &unitType)
        : name(elemName), n1(node1), n2(node2), value(val), unit(unitType) {

    }

    virtual string getInfo() const {
        ostringstream oss;
        oss << name << " " << n1->getName() << " " << n2->getName() << " " << value << " " << unit;
        return oss.str();
    }

    string getName () {
        return name;
    }
//    string getType () {
//        return type;
//    }

    virtual string getType() const = 0;

    virtual ~Element() {}
};

class Resistor : public Element {
    double conductance;
public:
    Resistor(const string &elemName, Node* node1, Node* node2, double resistance)
        : Element(elemName, node1, node2, resistance, "Ohm")
    {
        conductance = (resistance != 0) ? (1.0 / resistance) : 0;
    }

    string getType() const override { return "Resistor"; }

};

class Capacitor : public Element {
public:
    Capacitor(const string &elemName, Node* node1, Node* node2, double capacitance)
        : Element(elemName, node1, node2, capacitance, "Farad") {}

    string getType() const override { return "Capacitor"; }

};

class Inductor : public Element {
public:
    Inductor(const string &elemName, Node* node1, Node* node2, double inductance)
        : Element(elemName, node1, node2, inductance, "Henry") {}

    string getType() const override { return "Inductor"; }

};

class VoltageSource : public Element {
public:
    VoltageSource(const string &elemName, Node* node1, Node* node2, double voltage)
        : Element(elemName, node1, node2, voltage, "Volt") {
    };

    string getType() const override { return "VS"; }

};

class CurrentSource : public Element {
public:
    CurrentSource(const string &elemName, Node* node1, Node* node2, double current)
        : Element(elemName, node1, node2, current, "Ampere") {}

    string getType() const override { return "CS"; }

};
//-----------------------------------------------

class Circuit {
    vector<Node*> nodes;
    vector<Element*> elements;

public:
    ~Circuit() {
        for (auto node : nodes) delete node;
        for (auto elem : elements) delete elem;
    }

    Node* getOrCreateNode(const string& nodeName) {
        for (auto node : nodes) {
            if (node->getName() == nodeName)
                return node;
        }
        Node* newNode = new Node();
        newNode->getName() = nodeName;
        nodes.push_back(newNode);
        return newNode;
    }

    void addElement(Element* elem) {
        elements.push_back(elem);
    }

    void deleteElement(const string& elemName) {
        for (auto it = elements.begin(); it != elements.end(); ++it) {
            if ((*it)->getName() == elemName) {
                delete *it;
                elements.erase(it);
                break;
            }
        }
    }

    void listElements(const string& type = "") const {
        for (const auto& e : elements) {
            if (type == "" || e->getType() == type) {
                cout << e->getInfo() << endl;
            }
        }
    }

    void renameNode(const string& oldName, const string& newName) {
        Node* target = nullptr;
        for (auto node : nodes) {
            if (node->getName() == oldName)
                target = node;
            if (node->getName() == newName) {
                cout << "ERROR: Node name " << newName << " already exists\n";
                return;
            }
        }
        if (target) {
            target->getName() = newName;
            cout << "SUCCESS: Node renamed from " << oldName << " to " << newName << endl;
        } else {
            cout << "ERROR: Node " << oldName << " does not exist\n";
        }
    }

    void listNodes() const {
        cout << "Available nodes:\n";
        for (const auto& node : nodes) {
            cout << node->getName() << endl;
        }
    }
};






//--------------------- parse commands ---------------------
void processCommand(const string& cmd, Circuit& circuit) {
    smatch match;

    // ---------------------- Add Resistor ----------------------
    if (regex_match(cmd, match, regex(R"(add\s+R([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        string name = "R" + match[1].str();
        Node* n1 = circuit.getOrCreateNode(match[2].str());
        Node* n2 = circuit.getOrCreateNode(match[3].str());
        double val = parseNumber(match[4].str());
        circuit.addElement(new Resistor(name, n1, n2, val));
        return;
    }

    // ---------------------- Delete Resistor ----------------------
    if (regex_match(cmd, match, regex(R"(delete\s+R([^ ]+))"))) {
        string name = "R" + match[1].str();
        circuit.deleteElement(name);
        return;
    }

    // ---------------------- Add Capacitor ----------------------
    if (regex_match(cmd, match, regex(R"(add\s+C([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        string name = "C" + match[1].str();
        Node* n1 = circuit.getOrCreateNode(match[2].str());
        Node* n2 = circuit.getOrCreateNode(match[3].str());
        double val = parseNumber(match[4].str());
        circuit.addElement(new Capacitor(name, n1, n2, val));
        return;
    }

    // ---------------------- List Elements ----------------------
    if (regex_match(cmd, match, regex(R"(list\s+([A-Za-z]+))"))) {
        circuit.listElements(match[1].str());
        return;
    }

    if (regex_match(cmd, regex(R"(list)"))) {
        circuit.listElements();
        return;
    }

    // ---------------------- Rename Node ----------------------
    if (regex_match(cmd, match, regex(R"(rename\s+node\s+([^ ]+)\s+([^ ]+))"))) {
        circuit.renameNode(match[1].str(), match[2].str());
        return;
    }

    // ---------------------- Show Node List ----------------------
    if (regex_match(cmd, regex(R"(\.nodes)"))) {
        circuit.listNodes();
        return;
    }

    cout << "ERROR: Unknown or malformed command" << endl;
}






int main() {

    Circuit circuit;
    string cmd;

    cout << "Enter command (or 'exit'): \n";
    while (getline(cin, cmd)) {
        if (cmd == "exit") break;
        processCommand(cmd, circuit);
    }

    return 0;
}
