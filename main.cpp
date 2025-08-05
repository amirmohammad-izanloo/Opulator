#include <bits/stdc++.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL2_gfx.h>
#include <fstream>
#include <sstream>

using namespace std;

//================ Custom Exception Classes =================
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

class DuplicateGroundNodeException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Duplicate Ground Added.";
    }
};

class SingularMatrixException : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Zero pivot encountered. The system is singular.";
    }
};

//================ Helper Function ===========================
// Parses a number string with support for unit prefixes (k, m, u)
double parseNumber(string input) {
    double result = 0;
    char last = input[input.length() - 1];
    string numberPart;

    if (last == 'M') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) * 1000000; // Mega = 10^6
    } else if (last == 'n') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) / 1000000000; // nano = 10^-9
    } else if (last == 'k' || last == 'K') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) * 1000;
    } else if (last == 'm') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) / 1000;
    } else if (input.find('e') != string::npos) {
        result = stod(input); // Handle scientific notation
    } else if (last == 'u' || last == 'U') {
        numberPart = input.substr(0, input.length() - 1);
        result = stod(numberPart) / 1000000;
    } else {
        result = stod(input);
    }
    return result;
}

//================ Schematic Structure & Global Storage ============
struct Schematic {
    string name;
    string schematichPath;
    vector<string> lines;
};

vector<Schematic> gSchematics;

// Extracts a file name from a given file path
string extractFileName(const string &filePath) {
    size_t pos = filePath.find_last_of("/\\");
    if (pos != string::npos)
        return filePath.substr(pos+1);
    else
        return filePath;
}

//================ Gaussian Elimination Function =================

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        double pivot = A[i][i];
        if (abs(pivot) < 1e-12)
            throw SingularMatrixException();

        for (int j = i; j < n; j++)
            A[i][j] /= pivot;
        b[i] /= pivot;

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j < n; j++)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
    }
    return x;
}

//================ Circuit Classes (Nodes/Elements) ==============

class Node {
    string name;
    int x, y;
public:
    Node(const string& nodeName, int xCoord = 0, int yCoord = 0)
            : name(nodeName), x(xCoord), y(yCoord) {}

    string getName() const { return name; }
    void setName(const string& newName) { name = newName; }

    int getX() const { return x; }
    int getY() const { return y; }
    void setX(int xCoord) { x = xCoord; }
    void setY(int yCoord) { y = yCoord; }
};

class Element {
protected:
    string name;
    Node* n1;
    Node* n2;
    double value;
    string unit;
public:
    Element(const string &elemName, Node* node1, Node* node2, double val, const string &unitType)
            : name(elemName), n1(node1), n2(node2), value(val), unit(unitType) { }
    virtual ~Element() { }

    virtual string getType() const = 0;

    virtual string getInfo() const {
        ostringstream oss;
        oss << name << " (" << getType() << ") "
            << n1->getName() << " " << n2->getName() << " "
            << value << " " << unit;
        return oss.str();
    }

    string getName() const { return name; }
    double getValue() const { return value; }
    Node* getNode1() const { return n1; }
    Node* getNode2() const { return n2; }
};

class Resistor : public Element {
public:
    Resistor(const string &elemName, Node* a, Node* b, double resistance)
            : Element(elemName, a, b, resistance, "Ohm") {
        if (resistance <= 0)
            throw InvalidValueException();
    }
    string getType() const override { return "Resistor"; }
};

class Capacitor : public Element {
public:
    Capacitor(const string &elemName, Node* a, Node* b, double capacitance)
            : Element(elemName, a, b, capacitance, "Farad") {
        if (capacitance <= 0)
            throw InvalidValueException();
    }
    string getType() const override { return "Capacitor"; }
};

class Inductor : public Element {
public:
    Inductor(const string &elemName, Node* a, Node* b, double inductance)
            : Element(elemName, a, b, inductance, "Henry") {
        if (inductance <= 0)
            throw InvalidValueException();
    }
    string getType() const override { return "Inductor"; }
};

class VoltageSource : public Element {
public:

    double amplitude = 0;
    double frequency = 0;
    bool isSine = false;

    VoltageSource(const string &elemName, Node* node1, Node* node2, double offset,
                  double amp = 0, double freq = 0)
            : Element(elemName, node1, node2, offset, "Volt"),
              amplitude(amp), frequency(freq)
    {
        if (amp > 0 && freq > 0)
            isSine = true;
    }

    string getType() const override { return "VoltageSource"; }
};

class CurrentSource : public Element {
public:
    CurrentSource(const string &elemName, Node* a, Node* b, double current)
            : Element(elemName, a, b, current, "Ampere") { }
    string getType() const override { return "CurrentSource"; }
};

//================ Circuit Class ============================
class Circuit {
public:
    vector<Node*> nodes;
    vector<Element*> elements;
    string groundName = "";
    double timeStep = 0.001;
    int totalSteps = 1000;

    string schematicPath = "";

    ~Circuit() {
        for (auto node : nodes)
            delete node;
        for (auto elem : elements)
            delete elem;
    }

    void reset() {
        for (auto node : nodes)
            delete node;
        nodes.clear();

        for (auto elem : elements)
            delete elem;
        elements.clear();

        groundName = "";
        schematicPath = "";
    }

    Node* getOrCreateNode(const string& nodeName) {
        for (auto node : nodes)
            if (node->getName() == nodeName)
                return node;
        Node* newNode = new Node(nodeName);
        nodes.push_back(newNode);
        return newNode;
    }

    void addElement(Element* elem) {
        elements.push_back(elem);
    }

    void deleteElement(const string &elemName) {
        for (auto it = elements.begin(); it != elements.end(); it++) {
            if ((*it)->getName() == elemName) {
                delete *it;
                elements.erase(it);
                cout << "Deleted element: " << elemName << endl;
                return;
            }
        }
        cout << "Element " << elemName << " not found." << endl;
    }

    void listElements(const string &filter = "") const {
        cout << "\nElements in Circuit:" << endl;
        for (auto elem : elements) {
            if (filter == "" || elem->getType() == filter)
                cout << elem->getInfo() << endl;
        }
    }

    void renameNode(const string &oldName, const string &newName) {
        for (auto node : nodes) {
            if (node->getName() == newName)
                throw DuplicateNameException();
        }
        Node* target = nullptr;
        for (auto node : nodes)
            if (node->getName() == oldName)
                target = node;
        if (!target)
            throw runtime_error(string("Error: Node ") + oldName + " does not exist.");
        target->setName(newName);
        cout << "SUCCESS: Node renamed from " << oldName << " to " << newName << endl;
    }

    void listNodes() const {
        cout << "\nNodes in Circuit:" << endl;
        for (auto node : nodes)
            cout << node->getName() << endl;
    }

    //------------------ Analyze ------------------

    void setTransientParams(double dt, int steps) {
        if (dt <= 0 || steps <= 0) throw InvalidValueException();
        timeStep = dt;
        totalSteps = steps;
    }

   void solveNodalAnalysis() {
        bool foundGround = false;
        for (auto node : nodes) {
            if (node->getName() == groundName) { foundGround = true; break; }
        }
        if (!foundGround)
            throw NoGroundException();

        vector<Node*> unknownNodes;
        for (auto node : nodes)
            if (node->getName() != groundName)
                unknownNodes.push_back(node);
        int n = unknownNodes.size();
        if (n == 0) {
            cout << "All nodes are ground!" << endl;
            return;
        }

        map<string, int> nodeIndex;
        for (int i = 0; i < n; i++)
            nodeIndex[unknownNodes[i]->getName()] = i;
        vector<vector<double>> A(n, vector<double>(n, 0.0));
        vector<double> b(n, 0.0);

        for (auto elem : elements) {
            // Process Resistors
            Resistor* r = dynamic_cast<Resistor*>(elem);
            if (r != nullptr) {
                double g = 1.0 / r->getValue();
                string nodeA = r->getNode1()->getName();
                string nodeB = r->getNode2()->getName();
                if (nodeA != groundName && nodeB != groundName) {
                    int i = nodeIndex[nodeA], j = nodeIndex[nodeB];
                    A[i][i] += g;
                    A[j][j] += g;
                    A[i][j] -= g;
                    A[j][i] -= g;
                } else if (nodeA == groundName && nodeB != groundName) {
                    int j = nodeIndex[nodeB];
                    A[j][j] += g;
                } else if (nodeB == groundName && nodeA != groundName) {
                    int i = nodeIndex[nodeA];
                    A[i][i] += g;
                }
            }
            // Process Voltage Sources
            VoltageSource* vs = dynamic_cast<VoltageSource*>(elem);
            if (vs != nullptr) {
                double V = vs->getValue();
                const double G_big = 1e6;
                string nodeA = vs->getNode1()->getName();
                string nodeB = vs->getNode2()->getName();
                if (nodeA == groundName && nodeB != groundName) {
                    int j = nodeIndex[nodeB];
                    A[j][j] += G_big;
                    b[j] -= V * G_big;
                } else if (nodeB == groundName && nodeA != groundName) {
                    int i = nodeIndex[nodeA];
                    A[i][i] += G_big;
                    b[i] += V * G_big;
                } else {
                    cout << "Warning: Voltage source " << vs->getName()
                         << " is not connected to ground. Skipping it." << endl;
                }
            }
            // Process Inductors for DC analysis
            if (Inductor* ind = dynamic_cast<Inductor*>(elem)) {
                const double G_big = 1e6; // A very high conductance
                string nodeA = ind->getNode1()->getName();
                string nodeB = ind->getNode2()->getName();
                if (nodeA == groundName && nodeB != groundName) {
                    int j = nodeIndex[nodeB];
                    A[j][j] += G_big;
                } else if (nodeB == groundName && nodeA != groundName) {
                    int i = nodeIndex[nodeA];
                    A[i][i] += G_big;
                } else {
                    // If neither node is ground, then without modified nodal analysis this might cause a singular system.
                    cout << "Warning: Inductor " << ind->getName() << " not referenced to ground may lead to a singular matrix." << endl;
                }
                // Continue to next element.
                continue;
            }

            // Process Current Sources
            CurrentSource* cs = dynamic_cast<CurrentSource*>(elem);
            if (cs != nullptr) {
                double I = cs->getValue();
                string nodeA = cs->getNode1()->getName();
                string nodeB = cs->getNode2()->getName();
                if (nodeA != groundName && nodeB != groundName) {
                    b[nodeIndex[nodeA]] += I;
                    b[nodeIndex[nodeB]] -= I;
                } else if (nodeA == groundName && nodeB != groundName) {
                    b[nodeIndex[nodeB]] -= I;
                } else if (nodeB == groundName && nodeA != groundName) {
                    b[nodeIndex[nodeA]] += I;
                }
            }
        }

        vector<double> voltages = gaussianElimination(A, b);
        cout << "\nNodal Voltages (with ground node " << groundName << " = 0 V):" << endl;
        for (auto &p : nodeIndex)
            cout << p.first << " = " << voltages[p.second] << " V" << endl;
    }

    void simulateTransient() {
        bool foundGround = false;
        for (auto node : nodes)
            if (node->getName() == groundName) { foundGround = true; break; }
        if (!foundGround)
            throw NoGroundException();

        vector<Node*> unknownNodes;
        for (auto node : nodes)
            if (node->getName() != groundName)
                unknownNodes.push_back(node);
        int n = unknownNodes.size();
        if (n == 0) {
            cout << "No unknown nodes." << endl;
            return;
        }
        map<string, int> nodeIndex;
        for (int i = 0; i < n; ++i)
            nodeIndex[unknownNodes[i]->getName()] = i;
        map<string, double> prevVoltages;
        for (auto node : unknownNodes)
            prevVoltages[node->getName()] = 0.0;

        // Static map to hold inductor currents (companion model state)
        static map<string, double> inductorCurrents;

        cout << "--- Transient Simulation Start ---" << endl;
        for (int step = 0; step < totalSteps; ++step) {
            double time = step * timeStep;
            vector<vector<double>> A(n, vector<double>(n, 0.0));
            vector<double> b(n, 0.0);
            // Stamp each circuit element
            for (auto elem : elements) {
                string na = elem->getNode1()->getName();
                string nb = elem->getNode2()->getName();
                int i = (na != groundName) ? nodeIndex[na] : -1;
                int j = (nb != groundName) ? nodeIndex[nb] : -1;
                if (elem->getType() == "Resistor") {
                    double g = 1.0 / elem->getValue();
                    if (i != -1) A[i][i] += g;
                    if (j != -1) A[j][j] += g;
                    if (i != -1 && j != -1) {
                        A[i][j] -= g;
                        A[j][i] -= g;
                    }
                }
                else if (elem->getType() == "CurrentSource") {
                    double I = elem->getValue();
                    if (i != -1) b[i] += I;
                    if (j != -1) b[j] -= I;
                }
                else if (elem->getType() == "Capacitor") {
                    double C = elem->getValue();
                    double G = C / timeStep; // Equivalent conductance
                    double v_prev_a = (i != -1) ? prevVoltages[na] : 0.0;
                    double v_prev_b = (j != -1) ? prevVoltages[nb] : 0.0;
                    double Ieq = G * (v_prev_a - v_prev_b);
                    if (i != -1) {
                        A[i][i] += G;
                        b[i] += Ieq;
                    }
                    if (j != -1) {
                        A[j][j] += G;
                        b[j] -= Ieq;
                    }
                    if (i != -1 && j != -1) {
                        A[i][j] -= G;
                        A[j][i] -= G;
                    }
                }
                else if (elem->getType() == "Inductor") {
                    // Companion model for inductors using the backward Euler method.
                    double L = elem->getValue();
                    double G_ind = timeStep / L;
                    string indName = elem->getName();
                    if (inductorCurrents.find(indName) == inductorCurrents.end())
                        inductorCurrents[indName] = 0.0; // Initial current = 0
                    double I_prev = inductorCurrents[indName];
                    if (i != -1) {
                        A[i][i] += G_ind;
                        b[i] -= I_prev;
                    }
                    if (j != -1) {
                        A[j][j] += G_ind;
                        b[j] += I_prev;
                    }
                    if (i != -1 && j != -1) {
                        A[i][j] -= G_ind;
                        A[j][i] -= G_ind;
                    }
                }
                else if (elem->getType() == "VoltageSource") {
                    VoltageSource* vs = dynamic_cast<VoltageSource*>(elem);
                    if (!vs) continue;
                    double V = vs->getValue();
                    if (vs->isSine)
                        V = V + vs->amplitude * sin(2 * M_PI * vs->frequency * time);
                    const double G_big = 1e6;
                    if (na == groundName && j != -1) {
                        A[j][j] += G_big;
                        b[j] -= V * G_big;
                    } else if (nb == groundName && i != -1) {
                        A[i][i] += G_big;
                        b[i] += V * G_big;
                    } else if (i != -1 && j != -1) {
                        cout << "Skipping VS " << vs->getName() << " (not grounded)" << endl;
                    }
                }
            }

            vector<double> x = gaussianElimination(A, b);

            // Update inductor currents based on the newly computed node voltages.
            for (auto elem : elements) {
                if (elem->getType() == "Inductor") {
                    string na = elem->getNode1()->getName();
                    string nb = elem->getNode2()->getName();
                    int i = (na != groundName) ? nodeIndex[na] : -1;
                    int j = (nb != groundName) ? nodeIndex[nb] : -1;
                    double V_ind = 0.0;
                    if (i != -1) V_ind += x[i];
                    if (j != -1) V_ind -= x[j];
                    double G_ind = timeStep / elem->getValue();
                    string indName = elem->getName();
                    inductorCurrents[indName] = inductorCurrents[indName] + G_ind * V_ind;
                }
            }
            cout << "t = " << time << " s: ";
            for (int k = 0; k < n; ++k) {
                string nodeName = unknownNodes[k]->getName();
                double v = x[k];
                prevVoltages[nodeName] = v;
                cout << nodeName << " = " << v << " V  ";
            }
            cout << endl;
        }
        cout << "--- Transient Simulation Done ---" << endl;
    }

    void simulateTransientWithCurrents() {
        bool foundGround = false;
        for (auto node : nodes)
            if (node->getName() == groundName) { foundGround = true; break; }
        if (!foundGround)
            throw NoGroundException();

        vector<Node*> unknownNodes;
        for (auto node : nodes)
            if (node->getName() != groundName)
                unknownNodes.push_back(node);

        int n = unknownNodes.size();
        if (n == 0) {
            cout << "No unknown nodes." << endl;
            return;
        }

        map<string, int> nodeIndex;
        for (int i = 0; i < n; ++i)
            nodeIndex[unknownNodes[i]->getName()] = i;

        map<string, double> prevVoltages;
        for (auto node : unknownNodes)
            prevVoltages[node->getName()] = 0.0;

        static map<string, double> inductorCurrents;

        cout << "--- Transient Simulation with Currents ---" << endl;

        for (int step = 0; step < totalSteps; ++step) {
            double time = step * timeStep;
            vector<vector<double>> A(n, vector<double>(n, 0.0));
            vector<double> b(n, 0.0);

            for (auto elem : elements) {
                string na = elem->getNode1()->getName();
                string nb = elem->getNode2()->getName();
                int i = (na != groundName) ? nodeIndex[na] : -1;
                int j = (nb != groundName) ? nodeIndex[nb] : -1;

                if (elem->getType() == "Resistor") {
                    double g = 1.0 / elem->getValue();
                    if (i != -1) A[i][i] += g;
                    if (j != -1) A[j][j] += g;
                    if (i != -1 && j != -1) {
                        A[i][j] -= g;
                        A[j][i] -= g;
                    }
                } else if (elem->getType() == "CurrentSource") {
                    double I = elem->getValue();
                    if (i != -1) b[i] += I;
                    if (j != -1) b[j] -= I;
                } else if (elem->getType() == "Capacitor") {
                    double C = elem->getValue();
                    double G = C / timeStep;
                    double v_prev_a = (i != -1) ? prevVoltages[na] : 0.0;
                    double v_prev_b = (j != -1) ? prevVoltages[nb] : 0.0;
                    double Ieq = G * (v_prev_a - v_prev_b);

                    if (i != -1) {
                        A[i][i] += G;
                        b[i] += Ieq;
                    }
                    if (j != -1) {
                        A[j][j] += G;
                        b[j] -= Ieq;
                    }
                    if (i != -1 && j != -1) {
                        A[i][j] -= G;
                        A[j][i] -= G;
                    }
                } else if (elem->getType() == "Inductor") {
                    double L = elem->getValue();
                    double G_ind = timeStep / L;
                    string indName = elem->getName();
                    if (inductorCurrents.find(indName) == inductorCurrents.end())
                        inductorCurrents[indName] = 0.0;
                    double I_prev = inductorCurrents[indName];

                    if (i != -1) {
                        A[i][i] += G_ind;
                        b[i] -= I_prev;
                    }
                    if (j != -1) {
                        A[j][j] += G_ind;
                        b[j] += I_prev;
                    }
                    if (i != -1 && j != -1) {
                        A[i][j] -= G_ind;
                        A[j][i] -= G_ind;
                    }
                } else if (elem->getType() == "VoltageSource") {
                    VoltageSource* vs = dynamic_cast<VoltageSource*>(elem);
                    if (!vs) continue;
                    double V = vs->getValue();
                    if (vs->isSine)
                        V += vs->amplitude * sin(2 * M_PI * vs->frequency * time);
                    const double G_big = 1e6;
                    if (na == groundName && j != -1) {
                        A[j][j] += G_big;
                        b[j] -= V * G_big;
                    } else if (nb == groundName && i != -1) {
                        A[i][i] += G_big;
                        b[i] += V * G_big;
                    } else if (i != -1 && j != -1) {
                        A[i][i] += G_big;
                        A[j][j] += G_big;
                        A[i][j] -= G_big;
                        A[j][i] -= G_big;
                        b[i] += V * G_big;
                        b[j] -= V * G_big;
                    }
                }
            }

            vector<double> x = gaussianElimination(A, b);

            // به‌روزرسانی جریان سلف‌ها
            for (auto elem : elements) {
                if (elem->getType() == "Inductor") {
                    string na = elem->getNode1()->getName();
                    string nb = elem->getNode2()->getName();
                    int i = (na != groundName) ? nodeIndex[na] : -1;
                    int j = (nb != groundName) ? nodeIndex[nb] : -1;
                    double V_ind = 0.0;
                    if (i != -1) V_ind += x[i];
                    if (j != -1) V_ind -= x[j];
                    double G_ind = timeStep / elem->getValue();
                    string indName = elem->getName();
                    inductorCurrents[indName] += G_ind * V_ind;
                }
            }

            // ذخیره ولتاژها برای گام بعدی
            map<string, double> newVoltages;
            for (int k = 0; k < n; ++k) {
                newVoltages[unknownNodes[k]->getName()] = x[k];
            }

            cout << fixed << setprecision(6);
            cout << "t = " << time << " s:\n";

            for (int k = 0; k < n; ++k) {
                string nodeName = unknownNodes[k]->getName();
                cout << "  V(" << nodeName << ") = " << x[k] << " V\n";
            }

            for (auto elem : elements) {
                string na = elem->getNode1()->getName();
                string nb = elem->getNode2()->getName();
                double v_a = (na != groundName) ? newVoltages[na] : 0.0;
                double v_b = (nb != groundName) ? newVoltages[nb] : 0.0;
                double v_prev_a = (prevVoltages.count(na) ? prevVoltages[na] : 0.0);
                double v_prev_b = (prevVoltages.count(nb) ? prevVoltages[nb] : 0.0);
                double current = 0.0;

                if (elem->getType() == "Resistor") {
                    current = (v_a - v_b) / elem->getValue();
                } else if (elem->getType() == "Capacitor") {
                    double C = elem->getValue();
                    double dv = (v_a - v_b) - (v_prev_a - v_prev_b);
                    current = C * dv / timeStep;
                } else if (elem->getType() == "Inductor") {
                    current = inductorCurrents[elem->getName()];
                } else if (elem->getType() == "CurrentSource") {
                    current = elem->getValue();
                } else if (elem->getType() == "VoltageSource") {
                    current = NAN; // یا صرفاً چاپ نکنیم
                }

                if (!isnan(current))
                    cout << "  I(" << elem->getName() << ") = " << current << " A\n";
            }

            cout << endl;
            prevVoltages = newVoltages;
        }

        cout << "--- Transient Simulation with Currents Done ---" << endl;
    }

    // ------------------------- Ground ------------------------------

    void setGroundNode(const string& nodeName) {
        for (auto node : nodes) {
            if (node->getName() == nodeName) {
                groundName = nodeName;
                cout << "Ground node added: " << nodeName << endl;
                return;
            }
        }

        Node* newNode = new Node(nodeName);
        nodes.push_back(newNode);
        groundName = nodeName;
        cout << "Ground node created and added: " << nodeName << endl;
    }

    const string getGroundName() const {
        return groundName;
    }

};

Circuit circuit;
int menuLevel = 0;

void loadNewFile(const string &filePath) {
    ofstream infile(filePath, ios::app);
    if (!infile) {
        cout << "Error: Could not open file: " << filePath << endl;
        return;
    }

    bool exists = false;
    ifstream schList("C:\\\\Users\\\\Ared\\\\Desktop\\\\schList\\\\schList.txt");

    string line;
    while (getline(schList, line)) {
        if (line == filePath) {
            exists = true;
            break;
        }
    }
    schList.close();

    infile.close();
    if (!exists) {
        ofstream schematicsList("C:\\\\Users\\\\Ared\\\\Desktop\\\\schList\\\\schList.txt", ios::app);
        schematicsList << filePath << endl;
        schematicsList.close();
    }

    cout << "Schematic loaded: " << extractFileName(filePath) << endl;
}

void showSchematicsMenu() {
    gSchematics.clear();
    ifstream schematicsList("C:\\\\Users\\\\Ared\\\\Desktop\\\\schList\\\\schList.txt");
    string schematicPath;
    while(getline(schematicsList, schematicPath)) {
        ifstream infile(schematicPath);
        if (infile.is_open()) {
            Schematic sch;
            sch.name = extractFileName(schematicPath);
            sch.schematichPath = schematicPath;
            string line;
            while(getline(infile, line)) {
                // (Optionally: skip empty lines or comments)
                if (!line.empty())
                    sch.lines.push_back(line);
            }
            infile.close();
            gSchematics.push_back(sch);
        }
    }
    schematicsList.close();

    if(gSchematics.empty()) {
        cout << "No schematics available." << endl;
        return;
    }
    cout << "\nChoose existing schematic:" << endl;
    for (size_t i = 0; i < gSchematics.size(); i++) {
        cout << (i+1) << "- " << gSchematics[i].name << endl;
    }
    menuLevel = 1;
    cout << "Enter the schematic number (or type 'return' to go back):" << endl;
}

//================ Command Parsing Functions ====================
// Returns a vector of tokens; the first token identifies the command.
vector<string> parseCommandLine(const string &cmd) {
    vector<string> tokens;
    smatch match;

    if (menuLevel == 0 ) {

        // NewFile command: "NewFile <file_path>"
        if (regex_match(cmd, match, regex(R"(NewFile\s+([^ ]+))"))) {
            tokens.push_back("NewFile");
            tokens.push_back(match[1].str());
            return tokens;
        }
        // show existing schematics: exactly "-show existing schematics"
        if (regex_match(cmd, match, regex(R"(^-show\s+existing\s+schematics$)"))) {
            tokens.push_back("showSchematics");
            return tokens;
        }
    }

    if (menuLevel == 1) {
        // choose schematic: "<number>"
        if (regex_match(cmd, match, regex(R"(([^ ]+))"))) {
            tokens.push_back("chooseSchematic");
            tokens.push_back(match[1].str());
            return tokens;
        }
    }

    if (menuLevel == 2) {

        // exit existing schematics: exactly "-exit existing schematics"
        if (regex_match(cmd, match, regex(R"(^-exit\s+existing\s+schematics$)"))) {
            tokens.push_back("exitSchematics");
            return tokens;
        }
        // save existing schematics: exactly "-save existing schematics"
        if (regex_match(cmd, match, regex(R"(^-save\s+existing\s+schematics$)"))) {
            tokens.push_back("saveSchematics");
            return tokens;
        }
        // "return" command
        if (regex_match(cmd, regex(R"(return)"))) {
            tokens.push_back("return");
            return tokens;
        }
        // add resistor: "add R <Name> <node1> <node2> <value>"
        if (regex_match(cmd, match, regex(R"(add\s+R\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
            tokens.push_back("addR");
//            cout << "recivieng resistor" << endl;
            tokens.push_back(match[1].str());
            tokens.push_back(match[2].str());
            tokens.push_back(match[3].str());
            tokens.push_back(match[4].str());
            return tokens;
        }
        // add capacitor: "add C <Name> <node1> <node2> <value>"
        if (regex_match(cmd, match, regex(R"(add\s+C\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
            tokens.push_back("addC");
            tokens.push_back(match[1].str());
            tokens.push_back(match[2].str());
            tokens.push_back(match[3].str());
            tokens.push_back(match[4].str());
            return tokens;
        }
        // add inductor: "add L <Name> <node1> <node2> <value>"
        if (regex_match(cmd, match, regex(R"(add\s+L\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
            tokens.push_back("addL");
            tokens.push_back(match[1].str());
            tokens.push_back(match[2].str());
            tokens.push_back(match[3].str());
            tokens.push_back(match[4].str());
            return tokens;
        }
        // add voltage source: "add VS <Name> <node1> <node2> <value>"
        if (regex_match(cmd, match, regex(R"(add\s+VS\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
            tokens.push_back("addVS");
            tokens.push_back(match[1].str());
            tokens.push_back(match[2].str());
            tokens.push_back(match[3].str());
            tokens.push_back(match[4].str());
            return tokens;
        }
        // add current source: "add CS <Name> <node1> <node2> <value>"
        if (regex_match(cmd, match, regex(R"(add\s+CS\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+([^ ]+))"))) {
            tokens.push_back("addCS");
            tokens.push_back(match[1].str());
            tokens.push_back(match[2].str());
            tokens.push_back(match[3].str());
            tokens.push_back(match[4].str());
            return tokens;
        }
        // Voltage Source with SIN(...) syntax
        if (regex_match(cmd, match, regex(R"(add\s+V\s+([^ ]+)\s+([^ ]+)\s+([^ ]+)\s+SIN\(\s*([^ ]+)\s+([^ ]+)\s+([^ ]+)\s*\))"))) {
            tokens.push_back("addVS_SIN");
            for (int i = 1; i <= 6; ++i)
                tokens.push_back(match[i].str());
            return tokens;
        }
        if (regex_match(cmd, match, regex(R"(add\s+GND\s+([^ ]+))"))) {
            tokens.push_back("addGND");
            tokens.push_back(match[1].str()); // nodeName
            return tokens;
        }
        // delete and list commands are similar...
        if (regex_match(cmd, match, regex(R"(delete\s+([^ ]+))"))) {
            tokens.push_back("delete");
            tokens.push_back(match[1].str());
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
        if (regex_match(cmd, match, regex(R"(analyze)"))) {
            tokens.push_back("analyze");
            return tokens;
        }
        if (regex_match(cmd, match, regex(R"(\.step\s+([0-9eE\.\-]+)\s+(\d+))"))) {
            tokens.push_back(".step");
            tokens.push_back(match[1].str());  // dt
            tokens.push_back(match[2].str());  // number of steps
            return tokens;
        }
        if (regex_match(cmd, regex(R"(transient)"))) {
            tokens.push_back("transient");
            return tokens;
        }
        if (regex_match(cmd, regex(R"(transient_currents)"))) {
            tokens.push_back("transient_currents");
            return tokens;
        }
    }

    return tokens;
}

//================ Input Handler Function ====================

void chooseSchematic(const string &numStr);
void exitSchematic(Circuit &circuit);
void saveSchematic(Circuit &circuit);

void inputHandler(const string &input, Circuit &circuit) {
    vector<string> tokens = parseCommandLine(input);
    if (tokens.empty()) {
        cout << "ERROR: Unknown or malformed command" << endl;
        return;
    }
    string action = tokens[0];
    if (action == "NewFile") {
        string filePath = tokens[1];
        loadNewFile(filePath);
    }
    else if (action == "showSchematics") {
        showSchematicsMenu();
    }
    else if (action == "saveSchematics") {
        saveSchematic(circuit);
    }
    else if (action == "exitSchematics") {
        exitSchematic(circuit);
    }
    else if (action == "chooseSchematic") {
        if (tokens.size() < 2) {
            cout << "Error : Inappropriate input" << endl;
        } else {
            chooseSchematic(tokens[1]);
        }
    }
    else if (action == "return") {
        // Do nothing, simply return to main menu.
    }
    else if (action == "addR") {
        string name = tokens[1];

        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }

        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);

        circuit.addElement(new Resistor(name, n1, n2, val));
    }
    else if (action == "addC") {
        string name = tokens[1];
        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new Capacitor(name, n1, n2, val));
    }
    else if (action == "addL") {
        string name = tokens[1];
        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new Inductor(name, n1, n2, val));
    }
    else if (action == "addVS") {
        string name = tokens[1];
        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new VoltageSource(name, n1, n2, val));
    }
    else if (action == "addCS") {
        string name = tokens[1];
        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double val = parseNumber(tokens[4]);
        circuit.addElement(new CurrentSource(name, n1, n2, val));
    }
    else if (action == "addGND") {
        if (tokens.size() >= 2) {
            if (circuit.getGroundName() != "") {
                throw DuplicateGroundNodeException();
            }
            circuit.setGroundNode(tokens[1]);
        } else {
            cout << "ERROR: Missing node name for GND command" << endl;
        }
    }
    else if (action == "addVS_SIN") {
        // قالب: add V<name> <node+> <node-> SIN(<offset> <amplitude> <frequency>)
        string name = tokens[1];
        bool duplicate = false;
        for (auto elem : circuit.elements) {
            if (elem->getName() == name) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            throw DuplicateNameException();
        }
        Node* n1 = circuit.getOrCreateNode(tokens[2]);
        Node* n2 = circuit.getOrCreateNode(tokens[3]);
        double offset = parseNumber(tokens[4]);
        double amp = parseNumber(tokens[5]);
        double freq = parseNumber(tokens[6]);
        circuit.addElement(new VoltageSource(name, n1, n2, offset, amp, freq));
    }
    else if (action == "transient") {
        if (circuit.getGroundName() != "")
            circuit.simulateTransient(); // tokens[1] = ground node name
        else
            cout << "ERROR: Missing ground node in 'transient' command.\n";
    }
    else if (action == "transient_currents") {
        if (circuit.getGroundName() != "")
            circuit.simulateTransientWithCurrents();
        else
            cout << "ERROR: Missing ground node.\n";
    }
    else if (action == "delete") {
        string name = tokens[1];
        circuit.deleteElement(name);
    }
    else if (action == "list") {
        if (tokens.size() == 2)
            circuit.listElements(tokens[1]);
        else
            circuit.listElements();
    }
    else if (action == "rename") {
        if (tokens.size() >= 3)
            circuit.renameNode(tokens[1], tokens[2]);
    }
    else if (action == ".nodes") {
        circuit.listNodes();
    }
    else if (action == "analyze") {
        cout << circuit.groundName << endl;
            circuit.solveNodalAnalysis();
    }
    else if (action == ".step") {
        if (tokens.size() >= 3) {
            double dt = parseNumber(tokens[1]);
            int steps = stoi(tokens[2]);
            circuit.setTransientParams(dt, steps);
            cout << "Time step = " << dt << " s, total steps = " << steps << endl;
        }
    }
    else {
        cout << "ERROR: Unknown command action" << endl;
    }
}


void chooseSchematic(const string &numStr) {
    int choice = 0;
    try {
        choice = stoi(numStr);
    } catch (...) {
        cout << "Error : Inappropriate input" << endl;
        return;
    }
    if(choice < 1 || choice > static_cast<int>(gSchematics.size())) {
        cout << "Error : Inappropriate input" << endl;
        return;
    }

    // Display the chosen schematic's content.
    menuLevel = 2;
    Schematic sch = gSchematics[choice - 1];
    cout << "\nSchematic (" << sch.name << ") content:" << endl;
    for (const auto &line : sch.lines) {
        cout << line << endl;
    }

//    cout << sch.lines.size() << endl;
    for (const auto &line : sch.lines) {
        string line2 = "add " + line;
        inputHandler(line2, circuit);
    }

    circuit.schematicPath = sch.schematichPath;

}

void saveSchematic(Circuit &circuit) {
    vector<string> newLines;
    for (auto element : circuit.elements) {
        string newLine = "";
        bool isSin  = false;
        if (element->getType() == "Resistor") {
            newLine = "R ";
        } else if (element->getType() == "Capacitor") {
            newLine = "C ";

        } else if (element->getType() == "Inductor") {
            newLine = "L ";

        } else if (element->getType() == "VoltageSource") {
            VoltageSource* vs = dynamic_cast<VoltageSource*>(element);
            if (vs && vs->isSine) {
                isSin = true;
                ostringstream oss;
                oss << "V " << vs->getName() << " "
                    << vs->getNode1()->getName() << " "
                    << vs->getNode2()->getName() << " "
                    << "SIN(" << vs->getValue() << " " << vs->amplitude << " " << vs->frequency << ")";
                newLines.push_back(oss.str());
                continue;
            } else {
                newLine = "VS ";
            }
        } else if (element->getType() == "CurrentSource") {
            newLine = "CS ";
        }

        if(!isSin) {
            newLine = newLine + element->getName() + " " + element->getNode1()->getName() + " " + element->getNode2()->getName()
                      + " " + to_string(element->getValue());
        }

        newLines.push_back(newLine);
    }

    if (circuit.getGroundName() != "") {
        string newLine = "GND " + circuit.getGroundName();
        newLines.push_back(newLine);
    }

    ofstream outFile(circuit.schematicPath);
    if (!outFile) {
        cerr << "Schematic has been removed or transformed." << endl;
    } else {
        for (const auto& item : newLines) {
            outFile << item << std::endl;
        }
    }
    outFile.close();
    exitSchematic(circuit);

    cout << "Enter a new file or choose an old one." << endl;
}


void exitSchematic(Circuit &circuit) {
    circuit.reset();
    menuLevel = 0;
}


//================ Process Input Loop ========================
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


//================ Basic Drawing Functions ====================

enum ToolType { NONE, RESISTOR, CAPACITOR, INDUCTOR, VSOURCE, CSOURCE, GROUND, WIRE };
struct PlacedElement {
    ToolType type;
    int x, y;
    string name;
};
vector<PlacedElement> placedElements;
ToolType currentTool = NONE;


void drawConnector(SDL_Renderer* r, int x, int y) {
    SDL_Rect pin = {x - 2, y - 2, 5, 5};
    SDL_RenderFillRect(r, &pin);
}

void drawWireEnds(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    SDL_RenderDrawLine(r, x1 - 10, y1, x1, y1);
    SDL_RenderDrawLine(r, x2, y2, x2 + 10, y2);
    drawConnector(r, x1 - 10, y1);
    drawConnector(r, x2 + 10, y2);
}

void drawMenuBar(SDL_Renderer* renderer, TTF_Font* font) {
    SDL_Rect menuBar = {0, 0, 1280, 25};
    SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);
    SDL_RenderFillRect(renderer, &menuBar);

    SDL_Color textColor = {255, 255, 255};
    SDL_Surface* surfFile = TTF_RenderText_Solid(font, "File", textColor);
    SDL_Surface* surfEdit = TTF_RenderText_Solid(font, "Edit", textColor);
    SDL_Texture* texFile = SDL_CreateTextureFromSurface(renderer, surfFile);
    SDL_Texture* texEdit = SDL_CreateTextureFromSurface(renderer, surfEdit);

    SDL_Rect fileRect = {10, 4, surfFile->w, surfFile->h};
    SDL_Rect editRect = {70, 4, surfEdit->w, surfEdit->h};
    SDL_RenderCopy(renderer, texFile, nullptr, &fileRect);
    SDL_RenderCopy(renderer, texEdit, nullptr, &editRect);

    SDL_FreeSurface(surfFile);
    SDL_FreeSurface(surfEdit);
    SDL_DestroyTexture(texFile);
    SDL_DestroyTexture(texEdit);
}

void drawFileDropdown(SDL_Renderer* renderer, TTF_Font* font) {
    SDL_Rect bg = {10, 25, 100, 70};
    SDL_SetRenderDrawColor(renderer, 80, 80, 80, 255);
    SDL_RenderFillRect(renderer, &bg);

    SDL_Color white = {255, 255, 255};
    const char* options[] = {"New", "Open", "Save"};

    for (int i = 0; i < 3; ++i) {
        SDL_Surface* surf = TTF_RenderText_Solid(font, options[i], white);
        SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, surf);
        SDL_Rect dst = {15, 30 + i * 20, surf->w, surf->h};
        SDL_RenderCopy(renderer, tex, nullptr, &dst);
        SDL_FreeSurface(surf);
        SDL_DestroyTexture(tex);
    }
}

void drawResistor(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);
    const int segments = 6;
    SDL_Point points[segments + 2];
    double dx = (x2 - x1) / (double)(segments + 1);
    double dy = (y2 - y1) / (double)(segments + 1);
    points[0] = {x1, y1};
    for (int i = 1; i <= segments; ++i) {
        double offset = (i % 2 == 0) ? -5 : 5;
        points[i].x = x1 + i * dx;
        points[i].y = y1 + i * dy + offset;
    }
    points[segments + 1] = {x2, y2};
    SDL_RenderDrawLines(r, points, segments + 2);
}

void drawCapacitor(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);
    int mx = (x1 + x2) / 2, my = (y1 + y2) / 2;
    SDL_RenderDrawLine(r, x1, y1, mx - 5, my);
    SDL_RenderDrawLine(r, mx - 5, my - 10, mx - 5, my + 10);
    SDL_RenderDrawLine(r, mx + 5, my - 10, mx + 5, my + 10);
    SDL_RenderDrawLine(r, mx + 5, my, x2, y2);
}

void drawInductor(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);
    int mx = (x1 + x2) / 2;
    arcRGBA(r, mx - 10, y1, 5, 0, 180, 0, 0, 0, 255);
    arcRGBA(r, mx, y1, 5, 0, 180, 0, 0, 0, 255);
    arcRGBA(r, mx + 10, y1, 5, 0, 180, 0, 0, 0, 255);
    SDL_RenderDrawLine(r, x1, y1, mx - 15, y1);
    SDL_RenderDrawLine(r, mx + 15, y1, x2, y2);
}

void drawVoltageSource(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);
    int cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    circleRGBA(r, cx, cy, 15, 0, 0, 0, 255);
    SDL_RenderDrawLine(r, cx - 8, cy - 5, cx - 8, cy + 5);
    SDL_RenderDrawLine(r, cx - 12, cy, cx - 4, cy);
    SDL_RenderDrawLine(r, cx + 4, cy, cx + 12, cy);
}

void drawCurrentSource(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);
    int cx = (x1 + x2) / 2, cy = (y1 + y2) / 2;
    circleRGBA(r, cx, cy, 15, 0, 0, 0, 255);
    double a = atan2(y2 - y1, x2 - x1);
    int xt = cx + cos(a) * 10, yt = cy + sin(a) * 10;
    int xb1 = cx - cos(a) * 10 + sin(a) * 5;
    int yb1 = cy - sin(a) * 10 - cos(a) * 5;
    int xb2 = cx - cos(a) * 10 - sin(a) * 5;
    int yb2 = cy - sin(a) * 10 + cos(a) * 5;
    filledTrigonRGBA(r, xt, yt, xb1, yb1, xb2, yb2, 0, 0, 0, 255);
}

void drawGround(SDL_Renderer* r, int x, int y) {
    SDL_RenderDrawLine(r, x, y, x, y + 8);
    SDL_RenderDrawLine(r, x - 6, y + 8, x + 6, y + 8);
    SDL_RenderDrawLine(r, x - 4, y + 10, x + 4, y + 10);
    SDL_RenderDrawLine(r, x - 2, y + 12, x + 2, y + 12);
    drawConnector(r, x, y); // only one connector at top
}

void drawToolbar(SDL_Renderer* renderer) {
    SDL_Rect bar = {0, 25, 1280, 50};
    SDL_SetRenderDrawColor(renderer, 240, 240, 240, 255);
    SDL_RenderFillRect(renderer, &bar);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderDrawRect(renderer, &bar);

    int x = 30;
    drawResistor(renderer, x, 50, x + 30, 50); x += 60;
    drawCapacitor(renderer, x, 50, x + 30, 50); x += 60;
    drawInductor(renderer, x, 50, x + 30, 50); x += 60;
    drawVoltageSource(renderer, x, 50, x + 30, 50); x += 60;
    drawCurrentSource(renderer, x, 50, x + 30, 50); x += 60;
    drawGround(renderer, x + 15, 50); x += 60;
    SDL_RenderDrawLine(renderer, x + 10, 40, x + 10, 60); // یک خط عمودی به‌عنوان نماد سیم
    x += 30;
}

bool wireStartSelected = false;
int wireX1 = 0, wireY1 = 0;
vector<pair<SDL_Point, SDL_Point>> wires;

void drawWires(SDL_Renderer* r, const vector<pair<SDL_Point, SDL_Point>>& wires) {
    for (auto& w : wires) {
        aalineRGBA(r, w.first.x, w.first.y, w.second.x, w.second.y, 0, 0, 0, 255);
    }
}

void drawPlacedElements(SDL_Renderer* r) {
    for (auto& e : placedElements) {
        int x1 = e.x - 15, y1 = e.y;
        int x2 = e.x + 15, y2 = e.y;
        switch (e.type) {
            case RESISTOR: drawResistor(r, x1, y1, x2, y2); break;
            case CAPACITOR: drawCapacitor(r, x1, y1, x2, y2); break;
            case INDUCTOR: drawInductor(r, x1, y1, x2, y2); break;
            case VSOURCE: drawVoltageSource(r, x1, y1, x2, y2); break;
            case CSOURCE: drawCurrentSource(r, x1, y1, x2, y2); break;
            case GROUND: drawGround(r, e.x, e.y); break;
            default: break;
        }
    }
}

void drawPreviewElement(SDL_Renderer* r, int preX, int preY) {

    int previewX = preX, previewY = preY;

    if (currentTool == NONE) return;
    int x1 = previewX - 15, y1 = previewY;
    int x2 = previewX + 15, y2 = previewY;

    switch (currentTool) {
        case RESISTOR: drawResistor(r, x1, y1, x2, y2); break;
        case CAPACITOR: drawCapacitor(r, x1, y1, x2, y2); break;
        case INDUCTOR: drawInductor(r, x1, y1, x2, y2); break;
        case VSOURCE: drawVoltageSource(r, x1, y1, x2, y2); break;
        case CSOURCE: drawCurrentSource(r, x1, y1, x2, y2); break;
        case GROUND: drawGround(r, previewX, previewY); break;
        default: break;
    }
}


// ===================== Connectors ========================

struct Connector {
    int id;
    int x, y;
};

vector<Connector> allConnectors;
map<int, vector<int>> connectorGraph; // id → list of connected ids

map<int, string> connectorToNode;

int connectorIdCounter = 0;

map<int, string> connectorToNodeName;

int addConnector(SDL_Renderer* renderer, int x, int y) {
    for (const auto& c : allConnectors) {
        if (abs(c.x - x) <= 2 && abs(c.y - y) <= 2) return c.id;
    }
    allConnectors.push_back({connectorIdCounter, x, y});
    return connectorIdCounter++;
}

void registerConnectorConnection(int id1, int id2) {
    connectorGraph[id1].push_back(id2);
    connectorGraph[id2].push_back(id1);
}

void extractConnectorsFromPlacedElements(SDL_Renderer* renderer) {
    for (auto& e : placedElements) {
        int x1 = e.x - 25, x2 = e.x + 25;
        int y = e.y;
        if (e.type == GROUND) {
            int id = addConnector(renderer, e.x, e.y);
        } else {
            int id1 = addConnector(renderer, x1, y);
            int id2 = addConnector(renderer, x2, y);
        }
    }
}

void extractConnectorsFromWires(SDL_Renderer* renderer) {
    for (auto& w : wires) {
        int id1 = addConnector(renderer, w.first.x, w.first.y);
        int id2 = addConnector(renderer, w.second.x, w.second.y);
        registerConnectorConnection(id1, id2);
    }
}

void analyzeNodeConnections(SDL_Renderer* renderer) {
    extractConnectorsFromPlacedElements(renderer);
    extractConnectorsFromWires(renderer);

    for (int i = 0; i < allConnectors.size(); ++i) {
        filledCircleRGBA(renderer, allConnectors[i].x, allConnectors[i].y, 3, 0, 0, 255, 255);

    }

    map<int, int> parent;
    function<int(int)> find = [&](int u) {
        if (parent[u] != u) parent[u] = find(parent[u]);
        return parent[u];
    };
    auto unite = [&](int u, int v) {
        parent[find(u)] = find(v);
    };

    // initialize each connector as its own parent
    for (const auto& c : allConnectors) {
        parent[c.id] = c.id;
    }

    // unite connected connectors via wires
    for (const auto& entry : connectorGraph) {
        int u = entry.first;
        const vector<int>& neighbors = entry.second;
        for (int v : neighbors) {
            unite(u, v);
        }
    }

    // assign node names to connected components
    map<int, string> rootToNodeName;
    int nodeCounter = 1;
    for (const auto& c : allConnectors) {
        int root = find(c.id);
        if (!rootToNodeName.count(root)) {
            rootToNodeName[root] = "n" + to_string(nodeCounter++);
        }
        connectorToNode[c.id] = rootToNodeName[root];
    }
}

vector<string> generateAddCommands(SDL_Renderer* renderer) {
    vector<string> commands;
    set<string> addedNodes;

    for (const auto& c : allConnectors) {
        string nodeName = connectorToNode[c.id];
        if (addedNodes.insert(nodeName).second) {
            commands.push_back("add node " + nodeName + " " + to_string(c.x) + " " + to_string(c.y));
        }
    }

    map<ToolType, char> typeChar = {
            {RESISTOR, 'R'}, {CAPACITOR, 'C'}, {INDUCTOR, 'L'},
            {VSOURCE, 'V'}, {CSOURCE, 'I'}, {GROUND, 'G'}
    };

    map<ToolType, int> typeCounter;

    for (const auto& e : placedElements) {
        if (e.type == WIRE) continue;
        int x1 = e.x - 25, x2 = e.x + 25;
        int y = e.y;
        int id1 = (e.type == GROUND) ? addConnector(renderer, e.x, y) : addConnector(renderer, x1, y);
        int id2 = (e.type == GROUND) ? id1 : addConnector(renderer, x2, y);
        string node1 = connectorToNode[id1];
        string node2 = connectorToNode[id2];

        string name = typeChar[e.type] + to_string(++typeCounter[e.type]);
        string line = "add " + string(1, typeChar[e.type]) + " " + name + " " + node1;
        if (e.type != GROUND) line += " " + node2;

        // sample values:
        if (e.type == RESISTOR) line += " 1000";
        else if (e.type == CAPACITOR) line += " 1e-6";
        else if (e.type == INDUCTOR) line += " 0.01";
        else if (e.type == VSOURCE) line += " 5";
        else if (e.type == CSOURCE) line += " 0.1";

        commands.push_back(line);
    }

    return commands;
}

void saveCircuitToFile(SDL_Renderer* renderer, const string& filename) {
    analyzeNodeConnections(renderer);
    vector<string> cmds = generateAddCommands(renderer);
    ofstream out(filename);
    for (const auto& line : cmds) out << line << "\n";
    out.close();
}


//================ Main Function =============================
int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    TTF_Init();
    SDL_Window* window = SDL_CreateWindow("Circuit Visualizer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    TTF_Font* font = TTF_OpenFont("C:\\Windows\\Fonts\\Arial.ttf", 14);


    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
//    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    SDL_RenderPresent(renderer);

    bool showFileMenu = false;
    bool fileClicked = false;

    bool showPreview = false;
    int preX, preY;

    SDL_Event e;
    bool quit = false;


    while (!quit) {
        int mx, my;
        SDL_GetMouseState(&mx, &my);
        bool hoveringFile = (mx >= 10 && mx <= 60 && my <= 25);
        if (hoveringFile) fileClicked = true;
        bool hoveringDropdown = (mx >= 10 && mx <= 110 && my > 25 && my <= 95);



        if (!hoveringDropdown && !hoveringFile) {
            fileClicked = false;
        }

        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                quit = true;
            } else if (e.type == SDL_MOUSEMOTION) {
                showFileMenu = hoveringFile || (hoveringDropdown && fileClicked);

                if (currentTool != NONE && e.motion.y > 75) {
                    preX = e.motion.x;
                    preY = e.motion.y;
                    showPreview = true;
                } else {
                    showPreview = false;
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                cout << "x: " << e.button.x << " y: " << e.button.y << endl;
                if (hoveringDropdown && showFileMenu) {
                    int option = (my - 30) / 20;
                    switch(option) {
                        case 0:
                            // TODO: Implement New File
                            break;
                        case 1:
                            // TODO: Implement Open File
                            break;
                        case 2:
                            saveCircuitToFile(renderer, "C:\\\\Users\\\\Ared\\\\Desktop\\\\Circuits\\\\circuit1.txt");
                            break;
                    }
                    fileClicked = false;
                    showFileMenu = false;
                } else if (my >= 25 && my <= 75) {
                    if (mx >= 30 && mx <= 60) currentTool = RESISTOR;
                    else if (mx >= 90 && mx <= 120) currentTool = CAPACITOR;
                    else if (mx >= 150 && mx <= 180) currentTool = INDUCTOR;
                    else if (mx >= 210 && mx <= 240) currentTool = VSOURCE;
                    else if (mx >= 270 && mx <= 300) currentTool = CSOURCE;
                    else if (mx >= 330 && mx <= 360) currentTool = GROUND;
                    else if (mx >= 390 && mx <= 420) currentTool = WIRE;

                    if (mx >= 30 && mx <= 360) {
                        preX = 0; preY = 0;
                        showPreview = true;
                    }
                } else if (my > 75) {
                    if (currentTool == WIRE) {
                        if (!wireStartSelected) {
                            wireX1 = mx; wireY1 = my;
                            wireStartSelected = true;
                        } else {
                            wires.push_back({{wireX1, wireY1}, {mx, my}});
                            wireStartSelected = false;
                            currentTool = NONE;
                            showPreview = false;
                        }
                    } else if (currentTool != NONE) {
                        placedElements.push_back({currentTool, mx, my});
                        currentTool = NONE;
                        showPreview = false;
                    }
                }
            }
        }

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);
        drawMenuBar(renderer, font);
        drawToolbar(renderer);
        drawPlacedElements(renderer);
        if (showFileMenu) drawFileDropdown(renderer, font);
        if (showPreview) drawPreviewElement(renderer, preX, preY);

        analyzeNodeConnections(renderer);

        drawWires(renderer, wires);

        SDL_RenderPresent(renderer);
    }


    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    return 0;
}