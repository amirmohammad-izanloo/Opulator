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
    virtual Node* getNode1() const { return n1; }
    virtual Node* getNode2() const { return n2; }

    string getName () {
        return name;
    }

    virtual double getValue() {
        return value;
    }

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
//------------------- Gaussian Elimination Solver -------------------
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    // Forward elimination.
    for (int i = 0; i < n; i++) {
        double pivot = A[i][i];
        if (pivot == 0) {
            cerr << "Zero pivot encountered. System may be singular." << endl;
            exit(1);
        }
        for (int j = i; j < n; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    // Back substitution.
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
    return x;
}
//---------------------------------------------------
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
        Node* newNode = new Node(nodeName);
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
            // Since node name is private, call setName.
            target->setName(newName);
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

    // -----------------  Nodal Analysis Solver -----------------
    void solveNodalAnalysis(const string &groundName) {
    // Build a list of unknown nodes (all nodes except the ground).
    vector<Node*> unknownNodes;
    for (auto node : nodes) {  // nodes is assumed as a container of Node* in your Circuit class
        if (node->getName() != groundName)
            unknownNodes.push_back(node);
    }
    int n = unknownNodes.size();
    if(n == 0) {
        cout << "No unknown nodes to solve for (all nodes are ground?)" << endl;
        return;
    }
    // Map node name to matrix index.
    map<string, int> nodeIndex;
    for (int i = 0; i < n; i++) {
        nodeIndex[unknownNodes[i]->getName()] = i;
    }

    // Initialize matrix A (n x n) and vector b.
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    // Loop over all elements.
    for (auto elem : elements) {
        // Check for resistors.
        Resistor* r = dynamic_cast<Resistor*>(elem);
        if(r != nullptr) {
            double g = 1.0 / r->getValue(); // conductance
            string nodeA = r->getNode1()->getName();
            string nodeB = r->getNode2()->getName();
            if(nodeA != groundName && nodeB != groundName) {
                int i = nodeIndex[nodeA];
                int j = nodeIndex[nodeB];
                A[i][i] += g;
                A[j][j] += g;
                A[i][j] -= g;
                A[j][i] -= g;
            } else if(nodeA == groundName && nodeB != groundName) {
                int j = nodeIndex[nodeB];
                A[j][j] += g;
            } else if(nodeB == groundName && nodeA != groundName) {
                int i = nodeIndex[nodeA];
                A[i][i] += g;
            }
        }

        // Check for current sources.
        CurrentSource* cs = dynamic_cast<CurrentSource*>(elem);
        if(cs != nullptr) {
            double I = cs->getValue(); // positive: injection into node1, leaving node2.
            string nodeA = cs->getNode1()->getName();
            string nodeB = cs->getNode2()->getName();
            if(nodeA != groundName && nodeB != groundName) {
                b[nodeIndex[nodeA]] += I;
                b[nodeIndex[nodeB]] -= I;
            } else if(nodeA == groundName && nodeB != groundName) {
                b[nodeIndex[nodeB]] -= I;
            } else if(nodeB == groundName && nodeA != groundName) {
                b[nodeIndex[nodeA]] += I;
            }
        }

        // Check for capacitors: For DC analysis, treat as open-circuit (no effect).
        Capacitor* cap = dynamic_cast<Capacitor*>(elem);
        if(cap != nullptr) {
            // No contribution in DC (steady-state) analysis.
        }

        // Check for voltage sources.
        VoltageSource* vs = dynamic_cast<VoltageSource*>(elem);
        if(vs != nullptr) {
            // Simplified handling: if one of the nodes is ground, force the other node to the source voltage.
            // We add a very large conductance (G_big) into the diagonal entry. Otherwise, we'll just print a warning.
            double V = vs->getValue();
            const double G_big = 1e6; // A large conductance value.
            string nodeA = vs->getNode1()->getName();
            string nodeB = vs->getNode2()->getName();
            if(nodeA == groundName && nodeB != groundName) {
                int j = nodeIndex[nodeB];
                A[j][j] += G_big;
                b[j] += V * G_big;
            } else if(nodeB == groundName && nodeA != groundName) {
                int i = nodeIndex[nodeA];
                A[i][i] += G_big;
                b[i] += V * G_big;
            } else {
                // More advanced MNA needed if voltage source is between two non-ground nodes.
                cout << "Warning: Voltage source " << vs->getName()
                     << " requires full MNA. Skipping its contribution." << endl;
            }
        }

        // Check for inductors.
        Inductor* ind = dynamic_cast<Inductor*>(elem);
        if(ind != nullptr) {
            // In DC, an inductor behaves as a short circuit.
            // For nodal analysis, a perfect short means zero resistance.
            // However, to keep the matrix finite, we approximate it by a very high conductance.
            const double G_short = 1e6;
            string nodeA = ind->getNode1()->getName();
            string nodeB = ind->getNode2()->getName();
            if(nodeA != groundName && nodeB != groundName) {
                int i = nodeIndex[nodeA];
                int j = nodeIndex[nodeB];
                A[i][i] += G_short;
                A[j][j] += G_short;
                A[i][j] -= G_short;
                A[j][i] -= G_short;
            } else if(nodeA == groundName && nodeB != groundName) {
                int j = nodeIndex[nodeB];
                A[j][j] += G_short;
            } else if(nodeB == groundName && nodeA != groundName) {
                int i = nodeIndex[nodeA];
                A[i][i] += G_short;
            }
        }
    }

    // Solve A*x = b using Gaussian elimination.
    vector<double> voltages = gaussianElimination(A, b);

    cout << "\nNodal Voltages (with ground node " << groundName << " = 0 V):\n";
    for (auto &p : nodeIndex) {
        cout << p.first << " = " << voltages[p.second] << " V" << endl;
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
    // ---------------------- Analyze (Nodal Analysis) ----------------------
    if (regex_match(cmd, match, regex(R"(analyze\s+([^ ]+))"))) {
        string ground = match[1].str();
        circuit.solveNodalAnalysis(ground);
        return;
    }

    cout << "ERROR: Unknown or malformed command" << endl;
}





// -------------------- main function with a test -----------------------
int main() {
    Circuit circuit;
    // Pre-create nodes using getOrCreateNode.
    Node* n1 = circuit.getOrCreateNode("N1");
    Node* n2 = circuit.getOrCreateNode("N2");
    Node* n3 = circuit.getOrCreateNode("N3"); // ground

    // Add elements.
    circuit.addElement(new Resistor("R1", n1, n2, 1000));
    circuit.addElement(new Resistor("R2", n2, n3, 2000));
    circuit.addElement(new CurrentSource("I1", n1, n3, 0.01));

    circuit.listElements();
    circuit.listNodes();

    // Run nodal analysis with N3 as ground.
    circuit.solveNodalAnalysis("N3");

    return 0;
}
