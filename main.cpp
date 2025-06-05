#include <bits/stdc++.h>

using namespace std;

//------------------- Error exception ----------------------
class NoGroundException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: No ground node detected in the circuit.";
    }
};

class InvalidValueException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Negative or zero value for a component is invalid.";
    }
};

class DuplicateNameException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Duplicate node or component name detected.";
    }
};

class DependentSourceException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Dependent source has an undefined control element.";
    }
};

class MissingParameterException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Missing parameters for time-dependent source.";
    }
};

class SingularMatrixException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Zero pivot encountered. The system is singular.";
    }
};

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
    Element(const string &elemName, Node* node1, Node* node2, double val, const string &unitType)
            : name(elemName), n1(node1), n2(node2), value(val), unit(unitType) {}

    virtual string getInfo() const {
        ostringstream oss;
        oss << name << " " << n1->getName() << " " << n2->getName() << " " << value << " " << unit;
        return oss.str();
    }
    virtual Node* getNode1() const { return n1; }
    virtual Node* getNode2() const { return n2; }

    string getName() { return name; }

    virtual double getValue() { return value; }

    virtual string getType() const = 0;

    virtual ~Element() {}
};

class Resistor : public Element {
    double conductance;
public:
    Resistor(const string &elemName, Node* node1, Node* node2, double resistance)
            : Element(elemName, node1, node2, resistance, "Ohm") {
        if (resistance <= 0)
            throw InvalidValueException();
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
private:
    double current = 0.0; // جریان لحظه‌ای سلف

public:
    Inductor(const string &elemName, Node* node1, Node* node2, double inductance)
            : Element(elemName, node1, node2, inductance, "Henry") {}

    string getType() const override { return "Inductor"; }

    double getCurrent() const { return current; }
    void setCurrent(double newCurrent) { current = newCurrent; }
};


class VoltageSource : public Element {
public:
    VoltageSource(const string &elemName, Node* node1, Node* node2, double voltage)
            : Element(elemName, node1, node2, voltage, "Volt") {}

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
    for (int i = 0; i < n; i++) {
        double pivot = A[i][i];
        if (abs(pivot) < 1e-10) {
            throw SingularMatrixException();
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
    map<string, double> prevNodeVoltages;
    map<string, double> capacitorCurrents;
    map<string, double> inductorCurrents;
    set<string> groundNodes; // مجموعه گره‌های زمین
public:
    double timeStep = 0.001;
    double tStart = 0.0;
    double tStop = 0.1;
    double tMaxStep = 0.001;
    int totalSteps = 100;

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
        for (auto node : nodes) {
            if (node->getName() == newName)
                throw DuplicateNameException();
        }
        Node* target = nullptr;
        for (auto node : nodes) {
            if (node->getName() == oldName)
                target = node;
        }
        if (target) {
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

    void clearPrevVoltages() {
        prevNodeVoltages.clear();
        for (auto node : nodes) {
            prevNodeVoltages[node->getName()] = 0.0;
        }
    }

    void updateCapacitorCurrents() {
        capacitorCurrents.clear();
        for (auto elem : elements) {
            if (auto cap = dynamic_cast<Capacitor*>(elem)) {
                string name = cap->getName();
                double C = cap->getValue();
                double Vprev = prevNodeVoltages[cap->getNode1()->getName()] -
                               prevNodeVoltages[cap->getNode2()->getName()];
                double Geq = C / timeStep;
                capacitorCurrents[name] = -Geq * Vprev;
            }
        }
    }

    // grok #####################################################################
    void setGroundNode(const string& nodeName) {
        for (auto node : nodes) {
            if (node->getName() == nodeName) {
                groundNodes.insert(nodeName);
                cout << "Ground node added: " << nodeName << endl;
                return;
            }
        }
        // اگر گره وجود ندارد، آن را ایجاد می‌کنیم
        Node* newNode = new Node(nodeName);
        nodes.push_back(newNode);
        groundNodes.insert(nodeName);
        cout << "Ground node created and added: " << nodeName << endl;
    }

    const set<string>& getGroundNodes() const {
        return groundNodes;
    }

    void clearPrev() {
        prevNodeVoltages.clear();
        inductorCurrents.clear();
        for (Node* node : nodes) {
            prevNodeVoltages[node->getName()] = 0.0;
            for (auto elem : elements) {
                if (dynamic_cast<Inductor*>(elem)) {
                    inductorCurrents[elem->getName()] = 0.0;
                }
            }
        }
    }



    void updateInductorCurrents() {
        for (auto elem : elements) {
            if (auto ind = dynamic_cast<Inductor*>(elem)) {
                int n1 = ind->getNode1()->getVoltage();
                int n2 = ind->getNode2()->getVoltage();

                double voltage = n1 - n2;
                double L = ind->getValue(); // چون از Element ارث‌بری کرده، getValue همان مقدار سلف است
                double dt = timeStep;

                double di = (voltage / L) * dt;
                double newCurrent = ind->getCurrent() + di;
                ind->setCurrent(newCurrent);
            }
        }
    }







    void solveTransientAnalysis() {
        if (groundNodes.empty()) {
            groundNodes.insert("0"); // پیش‌فرض "0" اگر هیچ زمینی تعریف نشده
            cout << "No ground nodes defined, using default ground: 0" << endl;
        }

        bool foundGround = false;
        for (auto node : nodes) {
            if (groundNodes.count(node->getName()) > 0) {
                foundGround = true;
                break;
            }
        }
        if (!foundGround) throw NoGroundException();

        vector<Node*> unknownNodes;
        for (auto node : nodes) {
            if (groundNodes.count(node->getName()) == 0) // گره‌های غیرزمین
                unknownNodes.push_back(node);
        }
        int n = unknownNodes.size();
        if (n == 0) {
            cout << "No unknown nodes to solve for." << endl;
            return;
        }

        map<string, int> nodeIndex;
        for (int i = 0; i < n; i++) {
            nodeIndex[unknownNodes[i]->getName()] = i;
        }

        vector<vector<double>> results(totalSteps, vector<double>(n, 0.0));
        vector<double> times(totalSteps, 0.0);

        clearPrev();
        updateCapacitorCurrents();
        updateInductorCurrents();

        for (int step = 0; step < totalSteps; step++) {
            double t = tStart + step * timeStep;
            times[step] = t;

            vector<vector<double>> A(n, vector<double>(n, 0.0));
            vector<double> b(n, 0.0);

            for (auto elem : elements) {
                if (auto r = dynamic_cast<Resistor*>(elem)) {
                    double g = 1.0 / r->getValue();
                    string nodeA = r->getNode1()->getName();
                    string nodeB = r->getNode2()->getName();
                    if (groundNodes.count(nodeA) == 0 && groundNodes.count(nodeB) == 0) {
                        int i = nodeIndex[nodeA];
                        int j = nodeIndex[nodeB];
                        A[i][i] += g;
                        A[j][j] += g;
                        A[i][j] -= g;
                        A[j][i] -= g;
                    } else if (groundNodes.count(nodeA) > 0 && groundNodes.count(nodeB) == 0) {
                        int j = nodeIndex[nodeB];
                        A[j][j] += g;
                    } else if (groundNodes.count(nodeB) > 0 && groundNodes.count(nodeA) == 0) {
                        int i = nodeIndex[nodeA];
                        A[i][i] += g;
                    }
                }

                if (auto cs = dynamic_cast<CurrentSource*>(elem)) {
                    double I = cs->getValue();
                    string nodeA = cs->getNode1()->getName();
                    string nodeB = cs->getNode2()->getName();
                    if (groundNodes.count(nodeA) == 0 && groundNodes.count(nodeB) == 0) {
                        b[nodeIndex[nodeA]] += I;
                        b[nodeIndex[nodeB]] -= I;
                    } else if (groundNodes.count(nodeA) > 0 && groundNodes.count(nodeB) == 0) {
                        b[nodeIndex[nodeB]] -= I;
                    } else if (groundNodes.count(nodeB) > 0 && groundNodes.count(nodeA) == 0) {
                        b[nodeIndex[nodeA]] += I;
                    }
                }

                if (auto cap = dynamic_cast<Capacitor*>(elem)) {
                    double C = cap->getValue();
                    double Geq = C / timeStep;
                    double Ieq = capacitorCurrents[cap->getName()];
                    string nodeA = cap->getNode1()->getName();
                    string nodeB = cap->getNode2()->getName();
                    if (groundNodes.count(nodeA) == 0 && groundNodes.count(nodeB) == 0) {
                        int i = nodeIndex[nodeA];
                        int j = nodeIndex[nodeB];
                        A[i][i] += Geq;
                        A[j][j] += Geq;
                        A[i][j] -= Geq;
                        A[j][i] -= Geq;
                        b[i] += Ieq;
                        b[j] -= Ieq;
                    } else if (groundNodes.count(nodeA) > 0 && groundNodes.count(nodeB) == 0) {
                        int j = nodeIndex[nodeB];
                        A[j][j] += Geq;
                        b[j] -= Ieq;
                    } else if (groundNodes.count(nodeB) > 0 && groundNodes.count(nodeA) == 0) {
                        int i = nodeIndex[nodeA];
                        A[i][i] += Geq;
                        b[i] += Ieq;
                    }
                }

                if (auto vs = dynamic_cast<VoltageSource*>(elem)) {
                    double V = vs->getValue();
                    const double G_big = 1e6;
                    string nodeA = vs->getNode1()->getName();
                    string nodeB = vs->getNode2()->getName();
                    if (groundNodes.count(nodeA) > 0 && groundNodes.count(nodeB) == 0) {
                        int j = nodeIndex[nodeB];
                        A[j][j] += G_big;
                        b[j] += V * G_big;
                    } else if (groundNodes.count(nodeB) > 0 && groundNodes.count(nodeA) == 0) {
                        int i = nodeIndex[nodeA];
                        A[i][i] += G_big;
                        b[i] += V * G_big;
                    } else {
                        cout << "Warning: Voltage source " << vs->getName()
                             << " requires full MNA. Skipping." << endl;
                    }
                }

                if (auto ind = dynamic_cast<Inductor*>(elem)) {
                    double L = ind->getValue();
                    double Geq = timeStep / L;
                    double Ieq = inductorCurrents[ind->getName()];
                    string nodeA = ind->getNode1()->getName();
                    string nodeB = ind->getNode2()->getName();
                    if (groundNodes.count(nodeA) == 0 && groundNodes.count(nodeB) == 0) {
                        int i = nodeIndex[nodeA];
                        int j = nodeIndex[nodeB];
                        A[i][i] += Geq;
                        A[j][j] += Geq;
                        A[i][j] -= Geq;
                        A[j][i] -= Geq;
                        b[i] += Ieq;
                        b[j] -= Ieq;
                    } else if (groundNodes.count(nodeA) > 0 && groundNodes.count(nodeB) == 0) {
                        int j = nodeIndex[nodeB];
                        A[j][j] += Geq;
                        b[j] -= Ieq;
                    } else if (groundNodes.count(nodeB) > 0 && groundNodes.count(nodeA) == 0) {
                        int i = nodeIndex[nodeA];
                        A[i][i] += Geq;
                        b[i] += Ieq;
                    }
                }
            }

            try {
                vector<double> voltages = gaussianElimination(A, b);
                for (int i = 0; i < n; i++) {
                    results[step][i] = voltages[i];
                    prevNodeVoltages[unknownNodes[i]->getName()] = voltages[i];
                }
                updateCapacitorCurrents();
                for (auto elem : elements) {
                    if (auto ind = dynamic_cast<Inductor*>(elem)) {
                        string name = ind->getName();
                        double V = prevNodeVoltages[ind->getNode1()->getName()] -
                                   prevNodeVoltages[ind->getNode2()->getName()];
                        double Geq = timeStep / ind->getValue();
                        inductorCurrents[name] = Geq * V + inductorCurrents[name];
                    }
                }
            } catch (const SingularMatrixException &e) {
                cout << e.what() << " at time " << t << endl;
                return;
            }
        }

        cout << "\nTransient Analysis Results (Ground nodes: ";
        for (const auto& gnd : groundNodes) cout << gnd << " ";
        cout << "= 0 V):\n";
        cout << "Time(s)";
        for (auto node : unknownNodes) {
            cout << "\t" << node->getName() << "(V)";
        }
        cout << endl;
        for (int step = 0; step < totalSteps; step++) {
            cout << times[step];
            for (int i = 0; i < n; i++) {
                cout << "\t" << results[step][i];
            }
            cout << endl;
        }
    }

    void solveNodalAnalysis(const string &groundName) {
            if (groundNodes.empty()) {
                groundNodes.insert("0");
                cout << "No ground nodes defined, using default ground: 0" << endl;
            } else {

            // Build a list of unknown nodes (all nodes except the ground).
            vector<Node*> unknownNodes;
            bool foundGround = false;
            for (auto node : nodes) {  // nodes is your container of Node* in your Circuit class
                if (node->getName() == groundName) {
                    foundGround = true;
                    break;
                }
            }
            if (!foundGround) throw NoGroundException();
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
    }

};

//--------------------- parse commands ---------------------
vector<string> parseCommandLine(const string &cmd) {
    vector<string> tokens;
    smatch match;
    if (regex_match(cmd, match, regex(R"(add\s+R([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("addR");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(delete\s+R([^ ]+))"))) {
        tokens.push_back("deleteR");
        tokens.push_back(match[1].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(add\s+C([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("addC");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(add\s+L([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("addL");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(add\s+V([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("addV");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(add\s+I([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("addI");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(add\s+GND\s+([^ ]+))"))) {
        tokens.push_back("addGND");
        tokens.push_back(match[1].str()); // nodeName
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(list\s+([A-Za-z]+))"))) {
        tokens.push_back("list");
        tokens.push_back(match[1].str());
        return tokens;
    }
    if (regex_match(cmd, regex(R"(list)"))) {
        tokens.push_back("list");
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(rename\s+node\s+([^ ]+)\s+([^ ]+))"))) {
        tokens.push_back("rename");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        return tokens;
    }
    if (regex_match(cmd, regex(R"(\.nodes)"))) {
        tokens.push_back(".nodes");
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(analyze\s+([^ ]+))"))) {
        tokens.push_back("analyze");
        tokens.push_back(match[1].str());
        return tokens;
    }
    if (regex_match(cmd, match, regex(R"(TRAN\s+([^ ]+)\s+([^ ]+)\s*([^ ]*)\s*([^ ]*))"))) {
        tokens.push_back("TRAN");
        tokens.push_back(match[1].str());
        tokens.push_back(match[2].str());
        tokens.push_back(match[3].str());
        tokens.push_back(match[4].str());
        return tokens;
    }
    return tokens;
}

// ------------------- Input Handler Function -------------------
void inputHandler(const string &input, Circuit &circuit) {
    vector<string> tokens = parseCommandLine(input);
    if (tokens.empty()) {
        cout << "ERROR: Unknown or malformed command" << endl;
        return;
    }
    string action = tokens[0];
    if (action == "addR") {
        string name = "R" + tokens[1];
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new Resistor(name, n1, n2, val));
    } else if (action == "deleteR") {
        string name = "R" + tokens[1];
        circuit.deleteElement(name);
    } else if (action == "addC") {
        string name = "C" + tokens[1];
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new Capacitor(name, n1, n2, val));
    } else if (action == "addL") {
        string name = "L" + tokens[1];
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new Inductor(name, n1, n2, val));
    } else if (action == "addV") {
        string name = "V" + tokens[1];
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new VoltageSource(name, n1, n2, val));
    } else if (action == "addI") {
        string name = "I" + tokens[1];
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new CurrentSource(name, n1, n2, val));
    } else if (action == "addGND") {
        if (tokens.size() >= 2) {
            circuit.setGroundNode(tokens[1]);
        } else {
            cout << "ERROR: Missing node name for GND command" << endl;
        }
    } else if (action == "list") {
        if (tokens.size() == 2)
            circuit.listElements(tokens[1]);
        else
            circuit.listElements();
    } else if (action == "rename") {
        if (tokens.size() >= 3)
            circuit.renameNode(tokens[1], tokens[2]);
    } else if (action == ".nodes") {
        circuit.listNodes();
    } else if (action == "analyze") {
        if (tokens.size() >= 2)
            circuit.solveNodalAnalysis(tokens[1]);
    } else if (action == "TRAN") {
        if (tokens.size() >= 3) {
            circuit.timeStep = parseNumber(tokens[1]);
            circuit.tStop = parseNumber(tokens[2]);
            circuit.tStart = tokens.size() > 3 && !tokens[3].empty() ? parseNumber(tokens[3]) : 0.0;
            circuit.tMaxStep = tokens.size() > 4 && !tokens[4].empty() ? parseNumber(tokens[4]) : circuit.timeStep;
            circuit.totalSteps = static_cast<int>((circuit.tStop - circuit.tStart) / circuit.timeStep);
            circuit.solveTransientAnalysis();
        } else {
            cout << "ERROR: Insufficient parameters for TRAN command" << endl;
        }
    } else {
        cout << "ERROR: Unknown command action" << endl;
    }
}

void processInput(Circuit &circuit) {
    string input;
    cout << "Enter command (or 'exit' to quit):" << endl;
    while (getline(cin, input)) {
        if (input == "exit")
            break;
        try {
            inputHandler(input, circuit);
        }
        catch (const exception &e) {
            cout << e.what() << endl;
        }
    }
}

int main() {
    try {
        Circuit circuit;
//        testRCCircuit(circuit);
         processInput(circuit);
    }
    catch (const exception &e) {
        cout << "Unexpected exception: " << e.what() << endl;
    }
    return 0;
}