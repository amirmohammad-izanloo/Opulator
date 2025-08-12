#include <bits/stdc++.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL2_gfx.h>
#include <fstream>
#include <sstream>
#include <complex>


#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

#include "save_load.h"
#include "win_save_dialog.h"
#include "win_open_dialog.h"



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

// ===== Wave storage =====
struct WaveStore {
    vector<double> t;                        // time (s)
    map<string, vector<double>> V; // node voltages
    map<string, vector<double>> I; // element currents (optional)
    void clear(){ t.clear(); V.clear(); I.clear(); }
};

WaveStore gWaves;


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
        if (resistance < 0)
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

    bool computeNodalVoltages(map<string,double>& outV) const;

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
                double R = r->getValue();
                if (R <= 1e-9) R = 1e-9; // جلوگیری از تقسیم بر صفر
                double g = 1.0 / R;
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

    // قبلیِ اشتباه داخل کلاس:
// void simulateTransientCapture(Circuit& circuit, WaveStore& ws) {

    void simulateTransientCapture(WaveStore& ws) {
        ws.clear();

        // از اینجا به بعد هرجا گفتی circuit، بکنش this->
        bool foundGround = false;
        for (auto node : this->nodes)
            if (node->getName() == this->groundName) { foundGround = true; break; }
        if (!foundGround) throw NoGroundException();

        vector<Node*> unknownNodes;
        for (auto node : this->nodes)
            if (node->getName() != this->groundName)
                unknownNodes.push_back(node);
        int n = (int)unknownNodes.size();
        if (n == 0) return;

        map<string,int> nodeIndex;
        for (int i = 0; i < n; ++i) nodeIndex[unknownNodes[i]->getName()] = i;

        map<string,double> prevVoltages;
        for (auto node : unknownNodes) prevVoltages[node->getName()] = 0.0;

        static map<string,double> inductorCurrents;
        for (auto &kv : inductorCurrents) kv.second = 0.0;

        for (int step = 0; step < this->totalSteps; ++step) {
            double time = step * this->timeStep;

            vector<vector<double>> A(n, vector<double>(n, 0.0));
            vector<double> b(n, 0.0);

            for (auto elem : this->elements) {
                string na = elem->getNode1()->getName();
                string nb = elem->getNode2()->getName();
                int i = (na != this->groundName) ? nodeIndex[na] : -1;
                int j = (nb != this->groundName) ? nodeIndex[nb] : -1;

                if (elem->getType() == "Resistor") {
                    double g = 1.0 / elem->getValue();
                    if (i != -1) A[i][i] += g;
                    if (j != -1) A[j][j] += g;
                    if (i != -1 && j != -1) { A[i][j] -= g; A[j][i] -= g; }
                }
                else if (elem->getType() == "CurrentSource") {
                    double I = elem->getValue();
                    if (i != -1) b[i] += I;
                    if (j != -1) b[j] -= I;
                }
                else if (elem->getType() == "Capacitor") {
                    double C = elem->getValue();
                    double G = C / this->timeStep;
                    double v_prev_a = (i != -1) ? prevVoltages[na] : 0.0;
                    double v_prev_b = (j != -1) ? prevVoltages[nb] : 0.0;
                    double Ieq = G * (v_prev_a - v_prev_b);
                    if (i != -1) { A[i][i] += G; b[i] += Ieq; }
                    if (j != -1) { A[j][j] += G; b[j] -= Ieq; }
                    if (i != -1 && j != -1) { A[i][j] -= G; A[j][i] -= G; }
                }
                else if (elem->getType() == "Inductor") {
                    double L = elem->getValue();
                    double G_ind = this->timeStep / L;
                    string indName = elem->getName();
                    if (!inductorCurrents.count(indName)) inductorCurrents[indName] = 0.0;
                    double I_prev = inductorCurrents[indName];
                    if (i != -1) { A[i][i] += G_ind; b[i] -= I_prev; }
                    if (j != -1) { A[j][j] += G_ind; b[j] += I_prev; }
                    if (i != -1 && j != -1) { A[i][j] -= G_ind; A[j][i] -= G_ind; }
                }
                else if (elem->getType() == "VoltageSource") {
                    auto* vs = dynamic_cast<VoltageSource*>(elem);
                    if (!vs) continue;
                    double V = vs->getValue();
                    if (vs->isSine) V = V + vs->amplitude * sin(2*M_PI*vs->frequency*time);
                    const double G_big = 1e6;
                    if (na == this->groundName && j != -1) { A[j][j] += G_big; b[j] -= V*G_big; }
                    else if (nb == this->groundName && i != -1) { A[i][i] += G_big; b[i] += V*G_big; }
                }
            }

            vector<double> x = gaussianElimination(A, b);

            for (auto elem : this->elements) {
                if (elem->getType() == "Inductor") {
                    string na = elem->getNode1()->getName();
                    string nb = elem->getNode2()->getName();
                    int i = (na != this->groundName) ? nodeIndex[na] : -1;
                    int j = (nb != this->groundName) ? nodeIndex[nb] : -1;
                    double V_ind = 0.0;
                    if (i != -1) V_ind += x[i];
                    if (j != -1) V_ind -= x[j];
                    double G_ind = this->timeStep / elem->getValue();
                    string indName = elem->getName();
                    inductorCurrents[indName] = inductorCurrents[indName] + G_ind * V_ind;
                }
            }

            ws.t.push_back(time);
            for (int k = 0; k < n; ++k) {
                string nodeName = unknownNodes[k]->getName();
                ws.V[nodeName].push_back(x[k]);
                prevVoltages[nodeName] = x[k];
            }
            // در صورت نیاز ws.I[...] را هم اینجا پر کن
        }
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


    bool computeDCOP(map<string,double>& outV,
                     map<string,double>& outI) const;



};

// ======================= PATCH A — AC Sweep (fixed) =======================

inline string fmt6(double val) {
    ostringstream oss;
    oss << fixed << setprecision(6) << val;
    return oss.str();
}
// نوع عدد مختلط واحد
using cd = complex<double>;
static const cd j(0.0, 1.0);

// logspace دقیق: num>=2 نقطه بین start..stop (هر دو > 0)

// ===================== end of PATCH A — AC Sweep (fixed) ===================


// هر چیزی با R <= 1e-9 را «سیم» حساب کنیم
inline bool isWireResistor(const Element* e) {
    if (auto r = dynamic_cast<const Resistor*>(e)) {
        return r->getValue() <= 1e-9 + 1e-15;
    }
    return false;
}

// فرمت خروجی؛ پیش‌فرض 6 رقم بعد اعشار


string formatWithSI(double value, const char* unit = "") {
    if (!isfinite(value)) return "N/A";

    static const struct { const char* pfx; double mul; } prefixes[] = {
            {"G", 1e9}, {"M", 1e6}, {"k", 1e3}, {"", 1.0},
            {"m", 1e-3}, {"u", 1e-6}, {"n", 1e-9}, {"p", 1e-12}
    };

    double av = fabs(value);
    const char* pfx = "";
    double mul = 1.0;

    for (auto& t : prefixes) {
        if (av >= t.mul || (t.mul == 1e-12 && av < t.mul)) {
            pfx = t.pfx; mul = t.mul;
            break;
        }
    }

    double scaled = value / mul;

    ostringstream oss;
    oss << fixed << setprecision(6) << scaled;
    string s = oss.str();

    // حذف صفرهای اضافی انتها
    if (s.find('.') != string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }

    s += pfx;
    if (unit && unit[0]) s += unit;
    return s;
}

bool Circuit::computeDCOP(map<string,double>& outV,
                          map<string,double>& outI) const {
    // وجود گراند
    bool foundGround = false;
    for (auto node : nodes)
        if (node->getName() == groundName) { foundGround = true; break; }
    if (!foundGround) throw NoGroundException();

    // لیست گره‌های مجهول (به‌جز زمین)
    vector<Node*> unknown;
    for (auto node : nodes)
        if (node->getName() != groundName) unknown.push_back(node);

    int n = (int)unknown.size();
    outV.clear(); outI.clear();
    if (n == 0) return true;

    // ایندکس‌گذاری گره‌ها
    map<string,int> idx;
    for (int i=0;i<n;++i) idx[unknown[i]->getName()] = i;

    // ماتریس Ax=b هم‌سبکِ solveNodalAnalysis
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    for (auto elem : elements) {
        // Resistor
        if (auto* r = dynamic_cast<Resistor*>(elem)) {
            double R = r->getValue(); if (R <= 1e-9) R = 1e-9;
            double g = 1.0 / R;
            string a = r->getNode1()->getName();
            string c = r->getNode2()->getName();
            bool ia = (a != groundName), ic = (c != groundName);
            if (ia && ic) { int i = idx[a], j = idx[c]; A[i][i]+=g; A[j][j]+=g; A[i][j]-=g; A[j][i]-=g; }
            else if (!ia && ic) { int j = idx[c]; A[j][j]+=g; }
            else if (ia && !ic) { int i = idx[a]; A[i][i]+=g; }
            continue;
        }
        // Voltage source (روش هدایت بزرگ)
        if (auto* vs = dynamic_cast<VoltageSource*>(elem)) {
            double V = vs->getValue();               // DC offset
            const double G_big = 1e6;
            string a = vs->getNode1()->getName();
            string c = vs->getNode2()->getName();
            if (a == groundName && c != groundName) { int j = idx[c]; A[j][j]+=G_big; b[j] -= V*G_big; }
            else if (c == groundName && a != groundName) { int i = idx[a]; A[i][i]+=G_big; b[i] += V*G_big; }
            // در غیر این صورت فعلاً رها
            continue;
        }
        // Inductor (در DC: اتصال کوتاه معادل → هدایت بزرگ به زمین اگر یک سر به زمینه)
        if (auto* ind = dynamic_cast<Inductor*>(elem)) {
            const double G_big = 1e6;
            string a = ind->getNode1()->getName();
            string c = ind->getNode2()->getName();
            if (a == groundName && c != groundName) { int j = idx[c]; A[j][j]+=G_big; }
            else if (c == groundName && a != groundName) { int i = idx[a]; A[i][i]+=G_big; }
            continue;
        }
        // Current source
        if (auto* cs = dynamic_cast<CurrentSource*>(elem)) {
            double I = cs->getValue();
            string a = cs->getNode1()->getName();
            string c = cs->getNode2()->getName();
            if (a != groundName && c != groundName) { b[idx[a]] += I; b[idx[c]] -= I; }
            else if (a == groundName && c != groundName) { b[idx[c]] -= I; }
            else if (c == groundName && a != groundName) { b[idx[a]] += I; }
            continue;
        }
        // Capacitor در DC بازمدار است → اثری روی A ندارد
    }

    // حل
    vector<double> x = gaussianElimination(A, b);

    // ولتاژها
    for (auto& kv : idx) outV[kv.first] = x[kv.second];

    // جریان هر المان (تعریف مثبت از n1→n2)
    using std::numeric_limits;
    double NaN = numeric_limits<double>::quiet_NaN();

    for (auto elem : elements) {
        if (isWireResistor(elem)) {
            continue;
        }
        string name = elem->getName();
        string a = elem->getNode1()->getName();
        string c = elem->getNode2()->getName();
        double Va = (a == groundName) ? 0.0 : outV[a];
        double Vc = (c == groundName) ? 0.0 : outV[c];

        if (dynamic_cast<Resistor*>(elem)) {
            double R = elem->getValue(); if (R <= 1e-9) R = 1e-9;
            outI[name] = (Va - Vc) / R;
        } else if (dynamic_cast<CurrentSource*>(elem)) {
            outI[name] = elem->getValue();
        } else if (dynamic_cast<Capacitor*>(elem)) {
            outI[name] = 0.0;
        } else if (dynamic_cast<Inductor*>(elem)) {
            outI[name] = NaN; // DC ایده‌آل نامشخص
        } else if (dynamic_cast<VoltageSource*>(elem)) {
            outI[name] = NaN; // بدون MNA درنمیاد
        } else {
            outI[name] = NaN;
        }
    }

    return true;
}


bool Circuit::computeNodalVoltages(map<string,double>& outV) const {
    // همان منطق solveNodalAnalysis، ولی بدون چاپ و برگرداندن map
    bool foundGround = false;
    for (auto node : nodes) {
        if (node->getName() == groundName) { foundGround = true; break; }
    }
    if (!foundGround) throw NoGroundException();

    vector<Node*> unknownNodes;
    for (auto node : nodes)
        if (node->getName() != groundName)
            unknownNodes.push_back(node);

    int n = (int)unknownNodes.size();
    if (n == 0) { outV.clear(); return true; }

    map<string,int> nodeIndex;
    for (int i = 0; i < n; ++i)
        nodeIndex[unknownNodes[i]->getName()] = i;

    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    for (auto elem : elements) {
        // Resistor
        if (auto* r = dynamic_cast<Resistor*>(elem)) {
            double R = r->getValue();
            if (R <= 1e-9) R = 1e-9;
            double g = 1.0 / R;
            string a = r->getNode1()->getName();
            string bnode = r->getNode2()->getName();
            if (a != groundName && bnode != groundName) {
                int i = nodeIndex[a], j = nodeIndex[bnode];
                A[i][i] += g; A[j][j] += g; A[i][j] -= g; A[j][i] -= g;
            } else if (a == groundName && bnode != groundName) {
                int j = nodeIndex[bnode]; A[j][j] += g;
            } else if (bnode == groundName && a != groundName) {
                int i = nodeIndex[a]; A[i][i] += g;
            }
            continue;
        }
        // VoltageSource (روش هدایت بزرگ)
        if (auto* vs = dynamic_cast<VoltageSource*>(elem)) {
            double V = vs->getValue(); // برای DC: offset
            const double G_big = 1e6;
            string a = vs->getNode1()->getName();
            string bnode = vs->getNode2()->getName();
            if (a == groundName && bnode != groundName) {
                int j = nodeIndex[bnode]; A[j][j] += G_big; b[j] -= V * G_big;
            } else if (bnode == groundName && a != groundName) {
                int i = nodeIndex[a]; A[i][i] += G_big; b[i] += V * G_big;
            } else {
                // رفرنس زمین ندارد → از نظر DC ساده، صرف‌نظر (همان هشدار نسخهٔ چاپی)
            }
            continue;
        }
        // Inductor → در DC معادل اتصال کوتاه (با G_big)
        if (auto* ind = dynamic_cast<Inductor*>(elem)) {
            const double G_big = 1e6;
            string a = ind->getNode1()->getName();
            string bnode = ind->getNode2()->getName();
            if (a == groundName && bnode != groundName) {
                int j = nodeIndex[bnode]; A[j][j] += G_big;
            } else if (bnode == groundName && a != groundName) {
                int i = nodeIndex[a]; A[i][i] += G_big;
            } else {
                // ممکنه تکینگی بده؛ اینجا چیزی اضافه نمی‌کنیم (مثل هشدار قبلی)
            }
            continue;
        }
        // CurrentSource
        if (auto* cs = dynamic_cast<CurrentSource*>(elem)) {
            double I = cs->getValue();
            string a = cs->getNode1()->getName();
            string bnode = cs->getNode2()->getName();
            if (a != groundName && bnode != groundName) {
                b[nodeIndex[a]] += I;
                b[nodeIndex[bnode]] -= I;
            } else if (a == groundName && bnode != groundName) {
                b[nodeIndex[bnode]] -= I;
            } else if (bnode == groundName && a != groundName) {
                b[nodeIndex[a]] += I;
            }
            continue;
        }
    }

    vector<double> x = gaussianElimination(A, b);
    outV.clear();
    for (auto &kv : nodeIndex) outV[kv.first] = x[kv.second];
    return true;
}


Circuit circuit;


//================ Basic Drawing Functions ====================

enum ToolType { NONE, RESISTOR, CAPACITOR, INDUCTOR, VSOURCE, CSOURCE, GROUND, WIRE , VSIN  };
struct PlacedElement {
    ToolType type;
    int x, y;
    string name;
    double value;
    double amp = 0.0;     // دامنه ولتاژ (V)
    double freq = 0.0;    // فرکانس (Hz)
    double phase = 0.0;
};


// برای کنترل فوکوس بین فیلدها (۰=Name, 1=Amp, 2=Freq, 3=Phase)
int vacFocus = 0;
vector<PlacedElement> placedElements;
ToolType currentTool = NONE;

// --- Grid & Snap settings ---
// --- Grid & Snap settings ---
static const int GRID_SPACING     = 10;  // فاصله شبکه (px)
static const int GRID_DOT_SIZE    = 1;   // اندازه نقطه خاکستری
static const int SNAP_BOX_HALF    = 5;   // ناحیه جذب
static const int ELEMENT_HALF_LEN = 30;  // نیم‌طول بدنه (مضرب 10)
static const int PIN_STUB         = 10;  // طول استاب پین

const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB; // = 40




void drawConnector(SDL_Renderer* r, int x, int y) {
    SDL_Rect pin = {x - 2, y - 2, 5, 5};
    SDL_RenderFillRect(r, &pin);
}

void drawWireEnds(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    // x1 و x2 = سرهای بدنه
    SDL_RenderDrawLine(r, x1 - PIN_STUB, y1, x1, y1);
    SDL_RenderDrawLine(r, x2, y2, x2 + PIN_STUB, y2);

    SDL_Rect pin1 = { (x1 - PIN_STUB) - 2, y1 - 2, 5, 5 };
    SDL_Rect pin2 = { (x2 + PIN_STUB) - 2, y2 - 2, 5, 5 };
    SDL_RenderFillRect(r, &pin1);
    SDL_RenderFillRect(r, &pin2);
}


void drawMenuBar(SDL_Renderer* renderer, TTF_Font* font) {
    SDL_Rect menuBar = {0, 0, 1280, 25};
    SDL_SetRenderDrawColor(renderer, 100, 100, 100, 255);
    SDL_RenderFillRect(renderer, &menuBar);

    SDL_Color textColor = {255, 255, 255};
    SDL_Surface* surfFile = TTF_RenderText_Solid(font, "File", textColor);
    SDL_Surface* surfEdit = TTF_RenderText_Solid(font, "Edit", textColor);
    SDL_Surface* surfRun  = TTF_RenderText_Solid(font, "Run", textColor);

    SDL_Texture* texFile = SDL_CreateTextureFromSurface(renderer, surfFile);
    SDL_Texture* texEdit = SDL_CreateTextureFromSurface(renderer, surfEdit);
    SDL_Texture* texRun  = SDL_CreateTextureFromSurface(renderer, surfRun);

    SDL_Rect fileRect = {10, 4, surfFile->w, surfFile->h};
    SDL_Rect editRect = {70, 4, surfEdit->w, surfEdit->h};
    SDL_Rect runRect  = {130, 4, surfRun->w, surfRun->h};

    SDL_RenderCopy(renderer, texFile, nullptr, &fileRect);
    SDL_RenderCopy(renderer, texEdit, nullptr, &editRect);
    SDL_RenderCopy(renderer, texRun,  nullptr, &runRect);

    SDL_FreeSurface(surfFile);
    SDL_FreeSurface(surfEdit);
    SDL_FreeSurface(surfRun);
    SDL_DestroyTexture(texRun);
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
void drawACVoltageSource(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    // بدنه: از x1 تا x2 روی یک خط افقی
    // مرکز و شعاع دایره
    int cx = (x1 + x2) / 2;
    int cy = y1;
    int radius = 12;

    // پین‌ها و نقطه‌های اتصال مثل بقیه‌ی منابع
    drawWireEnds(r, x1, y1, x2, y2);

    // دایره
    SDL_SetRenderDrawColor(r, 0, 0, 0, 255);
    for (int w = 0; w < 360; w += 2) {
        double ra1 = (w)     * M_PI / 180.0;
        double ra2 = (w + 2) * M_PI / 180.0;
        int px1 = cx + (int)(radius * cos(ra1));
        int py1 = cy + (int)(radius * sin(ra1));
        int px2 = cx + (int)(radius * cos(ra2));
        int py2 = cy + (int)(radius * sin(ra2));
        SDL_RenderDrawLine(r, px1, py1, px2, py2);
    }

    // موج سینوسی کوچک داخل دایره
    int amp = radius / 3;           // دامنه‌ی کوچک
    int left = cx - radius + 3;
    int right = cx + radius - 3;
    int prevx = left;
    int prevy = cy + (int)(amp * sin((prevx - left) * 2.0 * M_PI / (right - left)));
    for (int x = left + 1; x <= right; ++x) {
        int y = cy + (int)(amp * sin((x - left) * 2.0 * M_PI / (right - left)));
        SDL_RenderDrawLine(r, prevx, prevy, x, y);
        prevx = x; prevy = y;
    }
}

void drawCurrentSource(SDL_Renderer* r, int x1, int y1, int x2, int y2) {
    drawWireEnds(r, x1, y1, x2, y2);

    int cx = (x1 + x2) / 2;
    int cy = (y1 + y2) / 2;
    circleRGBA(r, cx, cy, 15, 0, 0, 0, 255);

    // جهت جریان: از x1 به x2
    double angle = atan2(y2 - y1, x2 - x1);

    // طول فلش داخل دایره
    int lineLen = 10;
    int tipLen = 6;

    // محاسبه‌ی محل شروع و پایان پیکان (خط)
    int xStart = cx - cos(angle) * lineLen;
    int yStart = cy - sin(angle) * lineLen;
    int xTip   = cx + cos(angle) * lineLen;
    int yTip   = cy + sin(angle) * lineLen;

    // رسم خط فلش
    SDL_RenderDrawLine(r, xStart, yStart, xTip, yTip);

    // محاسبه‌ی مثلث نوک فلش
    int baseLeftX  = xTip - cos(angle) * tipLen - sin(angle) * tipLen / 2;
    int baseLeftY  = yTip - sin(angle) * tipLen + cos(angle) * tipLen / 2;

    int baseRightX = xTip - cos(angle) * tipLen + sin(angle) * tipLen / 2;
    int baseRightY = yTip - sin(angle) * tipLen - cos(angle) * tipLen / 2;

    // رسم مثلث نوک پیکان
    filledTrigonRGBA(r, xTip, yTip, baseLeftX, baseLeftY, baseRightX, baseRightY, 0, 0, 0, 255);
}

void drawGround(SDL_Renderer* r, int x, int y) {
    SDL_RenderDrawLine(r, x, y, x, y + 8);
    SDL_RenderDrawLine(r, x - 6, y + 8, x + 6, y + 8);
    SDL_RenderDrawLine(r, x - 4, y + 10, x + 4, y + 10);
    SDL_RenderDrawLine(r, x - 2, y + 12, x + 2, y + 12);
    drawConnector(r, x, y); // only one connector at top
}

void drawToolbar(SDL_Renderer* renderer, TTF_Font* font) {
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
    auto drawVACIcon = [&](int cx, int cy) {
        int x1 = cx - ELEMENT_HALF_LEN, y = cy;
        int x2 = cx + ELEMENT_HALF_LEN;
        drawACVoltageSource(renderer, x1, y, x2, y);
        SDL_Color black{0,0,0};
        SDL_Surface* s = TTF_RenderText_Solid(font, "VAC", black);
        if (s) {
            SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, s);
            SDL_Rect dst{ cx - s->w/2, cy + 18, s->w, s->h };
            SDL_RenderCopy(renderer, t, nullptr, &dst);
            SDL_DestroyTexture(t);
            SDL_FreeSurface(s);
        }
    };

    // فراخوانی در مختصات مناسب
    drawVACIcon(465, 50);

}

bool wireStartSelected = false;
int wireX1 = 0, wireY1 = 0;
vector<pair<SDL_Point, SDL_Point>> wires;

void drawWires(SDL_Renderer* r, const vector<pair<SDL_Point, SDL_Point>>& wires) {
    SDL_SetRenderDrawColor(r, 0, 0, 0, 255);
    for (auto& w : wires) {
        // draw wire line
        aalineRGBA(r, w.first.x, w.first.y, w.second.x, w.second.y, 0, 0, 0, 255);
        // draw 5x5 connector squares at each end
        SDL_Rect pin1 = { w.first.x - 2, w.first.y - 2, 5, 5 };
        SDL_Rect pin2 = { w.second.x - 2, w.second.y - 2, 5, 5 };
        SDL_RenderFillRect(r, &pin1);
        SDL_RenderFillRect(r, &pin2);
    }
}


void drawPlacedElements(SDL_Renderer* r) {
    for (auto& e : placedElements) {
        int x1 = e.x - ELEMENT_HALF_LEN, y1 = e.y;
        int x2 = e.x + ELEMENT_HALF_LEN, y2 = e.y;
        switch (e.type) {
            case RESISTOR:  drawResistor(r, x1, y1, x2, y2); break;
            case CAPACITOR: drawCapacitor(r, x1, y1, x2, y2); break;
            case INDUCTOR:  drawInductor(r, x1, y1, x2, y2); break;
            case VSOURCE:   drawVoltageSource(r, x1, y1, x2, y2); break;
            case CSOURCE:   drawCurrentSource(r, x1, y1, x2, y2); break;
            case GROUND:    drawGround(r, e.x, e.y); break;
            case VSIN:      drawACVoltageSource(r, x1, y1, x2, y2); break;
            default: break;
        }
    }
}


void drawPreviewElement(SDL_Renderer* r, int preX, int preY) {
    if (currentTool == NONE) return;
    int x1 = preX - ELEMENT_HALF_LEN, y1 = preY;
    int x2 = preX + ELEMENT_HALF_LEN, y2 = preY;

    switch (currentTool) {
        case RESISTOR:  drawResistor(r, x1, y1, x2, y2); break;
        case CAPACITOR: drawCapacitor(r, x1, y1, x2, y2); break;
        case INDUCTOR:  drawInductor(r, x1, y1, x2, y2); break;
        case VSOURCE:   drawVoltageSource(r, x1, y1, x2, y2); break;
        case CSOURCE:   drawCurrentSource(r, x1, y1, x2, y2); break;
        case GROUND:    drawGround(r, preX, preY); break;
        case VSIN:      drawACVoltageSource(r, x1, y1, x2, y2); break;
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



int addConnector(int x, int y) {
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

void extractConnectorsFromPlacedElements() {
    for (auto& e : placedElements) {
        if (e.type == GROUND) {
            addConnector(e.x, e.y);
        } else {
            // پین واقعی = بدنه ± (ELEMENT_HALF_LEN) به‌علاوه استاب پین
            int xL = e.x - PIN_OFFSET;
            int xR = e.x + PIN_OFFSET;
            int y  = e.y;
            addConnector(xL, y);
            addConnector(xR, y);
        }
    }
}



void extractConnectorsFromWires() {
    // For graphical wires only: register endpoints but do NOT electrically connect them
    for (auto& w : wires) {
        int id1 = addConnector(w.first.x, w.first.y);
        int id2 = addConnector(w.second.x, w.second.y);
        registerConnectorConnection(id1, id2);
        // Skipped registerConnectorConnection to keep wire endpoints as separate nodes
    }
}

std::unordered_map<int, std::string> loadedSeed;  // ✅

void analyzeNodeConnections(const std::unordered_map<int,std::string>* seed = nullptr) {
    // 0) استخراج کانکتورها از صحنه
    extractConnectorsFromPlacedElements();
    extractConnectorsFromWires();

    // 1) DSU (Union-Find) با path compression
    std::map<int,int> parent;
    auto find = [&](int u){
        int r = u;
        while (parent[r] != r) r = parent[r];
        // compression
        while (parent[u] != u) { int p = parent[u]; parent[u] = r; u = p; }
        return r;
    };
    auto unite = [&](int a, int b){
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[ra] = rb;
    };

    // 2) هر کانکتور، والد خودش
    for (const auto& c : allConnectors) parent[c.id] = c.id;

    // 3) اگر seed داریم: کانکتورهای هم‌نام را با هم یکی کن
    if (seed && !seed->empty()) {
        std::unordered_map<std::string, std::vector<int>> buckets;
        buckets.reserve(seed->size());
        for (const auto& kv : *seed) {
            // kv.first = connectorId, kv.second = nodeName
            buckets[kv.second].push_back(kv.first);
        }
        for (auto& b : buckets) {
            const auto& ids = b.second;
            for (size_t i = 1; i < ids.size(); ++i) unite(ids[0], ids[i]);
        }
    }

    // 4) یگان‌سازی براساس سیم‌ها (گراف اتصالات)
    for (const auto& entry : connectorGraph) {
        int u = entry.first;
        const std::vector<int>& neighbors = entry.second;
        for (int v : neighbors) unite(u, v);
    }

    // 5) نام‌گذاری مؤلفه‌ها
    std::map<int, std::string> rootToNodeName; // root -> name
    std::unordered_set<std::string> usedNames; usedNames.reserve(256);

    // 5.a اگر seed داریم، تا حد ممکن نام را از seed بردار
    if (seed && !seed->empty()) {
        // برای هر کانکتور seed، نام را روی ریشه‌اش ست کن (اگر چند نام روی یک ریشه افتاد، اولی می‌ماند)
        for (const auto& kv : *seed) {
            int root = find(kv.first);
            auto it = rootToNodeName.find(root);
            if (it == rootToNodeName.end()) {
                rootToNodeName[root] = kv.second;
                usedNames.insert(kv.second);
            }
        }
    }

    // 5.b برای مؤلفه‌های بی‌نام، نام یکتای n1,n2,… بساز (بدون تداخل با usedNames)
    auto genName = [&](int& counter){
        std::string s;
        do { s = "n" + std::to_string(counter++); } while (usedNames.count(s));
        usedNames.insert(s);
        return s;
    };
    int nodeCounter = 1;

    for (const auto& c : allConnectors) {
        int root = find(c.id);
        if (!rootToNodeName.count(root)) {
            rootToNodeName[root] = genName(nodeCounter);
        }
    }

    // 6) پر کردن connectorToNode برای همهٔ کانکتورها
    connectorToNode.clear();
    for (const auto& c : allConnectors) {
        int root = find(c.id);
        connectorToNode[c.id] = rootToNodeName[root];
    }
}


// ---------------- Canvas → Circuit ----------------
void buildCircuit(Circuit& circuit) {
    // 0) پاک‌سازی مدار و بازسازی توپولوژی از روی بوم
    circuit.reset();
//    allConnectors.clear();
    connectorToNode.clear();
    connectorGraph.clear();

    if (loadedSeed.empty()) {
        // حالت معمول از روی بوم
        allConnectors.clear();
        connectorIdCounter = 0;
        analyzeNodeConnections(/*seed=*/nullptr);
    } else {
        // بلافاصله بعد از لود از فایل
        // (allConnectors قبلاً در openAndLoadToCanvas ساخته شده!)
        analyzeNodeConnections(&loadedSeed);
    }


    // 1) ساخت گره‌ها (با مختصات، اگر خواستی نگه داری)
    unordered_map<string, pair<int,int>> nodePos;
    for (const auto& c : allConnectors) {
        auto it = connectorToNode.find(c.id);
        if (it == connectorToNode.end()) continue;
        const string& nn = it->second;
        if (!nodePos.count(nn)) nodePos[nn] = {c.x, c.y};
    }
    for (auto& kv : nodePos) {
        Node* n = circuit.getOrCreateNode(kv.first);
        n->setX(kv.second.first);
        n->setY(kv.second.second);
    }

    // 2) افزودن المان‌ها بر پایهٔ نوع ابزار
    const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB;
    map<ToolType,int> typeCounter;
    auto autoName = [&](const PlacedElement& e)->string{
        static map<ToolType,char> L = {
                {RESISTOR,'R'},{CAPACITOR,'C'},{INDUCTOR,'L'},
                {VSOURCE,'V'}, {CSOURCE,'I'},  {GROUND,'G'},
                {VSIN,'V'}
        };
        if (!e.name.empty()) return e.name;
        return string(1, L[e.type]) + to_string(++typeCounter[e.type]);
    };

    for (const auto& e : placedElements) {
        if (e.type == WIRE) continue;

        // id کانکتورها دقیقاً با همان هندسه‌ی رسم
        int id1, id2;
        int y = e.y;
        if (e.type == GROUND) {
            id1 = addConnector(e.x, y);
            id2 = id1;
        } else {
            id1 = addConnector(e.x - PIN_OFFSET, y);
            id2 = addConnector(e.x + PIN_OFFSET, y);
        }

        auto it1 = connectorToNode.find(id1);
        auto it2 = connectorToNode.find(id2);
        if (it1 == connectorToNode.end() || it2 == connectorToNode.end()) continue;

        const string node1 = it1->second;
        const string node2 = it2->second;

        Node* n1 = circuit.getOrCreateNode(node1);
        Node* n2 = circuit.getOrCreateNode(node2);
        string name = autoName(e);

        switch (e.type) {
            case RESISTOR:
                circuit.addElement(new Resistor(name, n1, n2, e.value));
                break;
            case CAPACITOR:
                circuit.addElement(new Capacitor(name, n1, n2, e.value));
                break;
            case INDUCTOR:
                circuit.addElement(new Inductor(name, n1, n2, e.value));
                break;
            case VSOURCE: // منبع ولتاژ DC
                circuit.addElement(new VoltageSource(name, n1, n2, e.value));
                break;
            case VSIN:    // VAC: offset, amp, freq
                circuit.addElement(new VoltageSource(name, n1, n2, e.value, e.amp, e.freq));
                break;
            case CSOURCE:
                circuit.addElement(new CurrentSource(name, n1, n2, e.value));
                break;
            case GROUND:
                if (circuit.getGroundName().empty())
                    circuit.setGroundNode(node1);
                break;
            default: break;
        }
    }

    // 3) (اختیاری) تبدیل Wire های گرافیکی به R=0 اگر دو نود متفاوت را پل بزنند
    int widx = 1;
    for (const auto& w : wires) {
        int id1 = addConnector(w.first.x,  w.first.y);
        int id2 = addConnector(w.second.x, w.second.y);
        auto it1 = connectorToNode.find(id1);
        auto it2 = connectorToNode.find(id2);
        if (it1 == connectorToNode.end() || it2 == connectorToNode.end()) continue;
        const string& a = it1->second;
        const string& b = it2->second;
        if (a == b) continue;
        Node* n1 = circuit.getOrCreateNode(a);
        Node* n2 = circuit.getOrCreateNode(b);
        circuit.addElement(new Resistor("W" + to_string(widx++), n1, n2, 0.0));
    }
}


vector<string> generateAddCommands() {
    vector<string> commands;
    set<string> addedNodes;

    // 0) dump unique nodes (so loader can restore coordinates)
    for (const auto& c : allConnectors) {
        auto it = connectorToNode.find(c.id);
        if (it == connectorToNode.end()) continue;
        const string& nodeName = it->second;
        if (addedNodes.insert(nodeName).second) {
            commands.push_back("node " + nodeName + " " + to_string(c.x) + " " + to_string(c.y));
        }
    }

    // shared geometry (must match draw/extract)
    const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB;

    // type -> letter
    map<ToolType, char> typeChar = {
            {RESISTOR,'R'}, {CAPACITOR,'C'}, {INDUCTOR,'L'},
            {VSOURCE,'V'},  {CSOURCE,'I'},  {GROUND,'G'},
            {VSIN,'V'}
    };
    map<ToolType,int> typeCounter;

    // 1) elements
    for (const auto& e : placedElements) {
        if (e.type == WIRE) continue; // wires handled below (optional)

        // resolve connector ids via your own registry (exact IDs)
        int y  = e.y;
        int id1, id2;
        if (e.type == GROUND) {
            id1 = addConnector(e.x, y);
            id2 = id1;
        } else {
            id1 = addConnector(e.x - PIN_OFFSET, y);
            id2 = addConnector(e.x + PIN_OFFSET, y);
        }

        auto it1 = connectorToNode.find(id1);
        auto it2 = connectorToNode.find(id2);
        if (it1 == connectorToNode.end() || it2 == connectorToNode.end()) continue; // safety

        const string node1 = it1->second;
        const string node2 = it2->second;

        // element name: prefer user-provided; else auto
        string name = e.name;
        if (name.empty()) name = string(1, typeChar[e.type]) + to_string(++typeCounter[e.type]);

        // build line
        if (e.type == GROUND) {
            commands.push_back("G " + node1);
            continue;
        }

        // VSIN: SIN(offset amp freq)
        if (e.type == VSIN) {
            ostringstream oss;
            // e.value = DC offset, e.amp = amplitude, e.freq = frequency
            oss << "V " << name << " " << node1 << " " << node2
                << " SIN(" << e.value << " " << e.amp << " " << e.freq << ")";
            commands.push_back(oss.str());
            continue;
        }

        // generic 2‑terminal (R/C/L, VSOURCE DC, CSOURCE DC)
        ostringstream line;
        line << typeChar[e.type] << " " << name << " " << node1 << " " << node2 << " " << e.value;

        // current source direction is node1 -> node2 (kept as-is)
        commands.push_back(line.str());
    }

    // 2) optional: export graphical wires as R=0 only if they bridge different nodes
    int widx = 1;
    for (const auto& w : wires) {
        int id1 = addConnector(w.first.x,  w.first.y);
        int id2 = addConnector(w.second.x, w.second.y);
        auto it1 = connectorToNode.find(id1);
        auto it2 = connectorToNode.find(id2);
        if (it1 == connectorToNode.end() || it2 == connectorToNode.end()) continue;

        const string& n1 = it1->second;
        const string& n2 = it2->second;
        if (n1 == n2) continue; // skip self-loop wires
        commands.push_back("R W" + to_string(widx++) + " " + n1 + " " + n2 + " 0");
    }

    return commands;
}

bool isGroundPlaced = false;


void showDCOPWindow(const map<string,double>& Vmap,
                    const map<string,double>& Imap,
                    const string& gndName) {
    SDL_Window* w = SDL_CreateWindow("DC Operating Point",
                                     SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                     760, 520, SDL_WINDOW_SHOWN);
    SDL_Renderer* r = SDL_CreateRenderer(w, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* f = TTF_OpenFont("C:\\Windows\\Fonts\\Arial.ttf", 14);

    auto drawText = [&](const string& s, int x, int y){
        if (!f) return;
        SDL_Color black{0,0,0};
        SDL_Surface* surf = TTF_RenderText_Solid(f, s.c_str(), black);
        if (!surf) return;
        SDL_Texture* tex = SDL_CreateTextureFromSurface(r, surf);
        SDL_Rect dst{ x, y, surf->w, surf->h };
        SDL_RenderCopy(r, tex, nullptr, &dst);
        SDL_DestroyTexture(tex); SDL_FreeSurface(surf);
    };

    // مرتب‌سازی برای نمایش منظم
    vector<pair<string,double>> nodes(Vmap.begin(), Vmap.end());
    sort(nodes.begin(), nodes.end(), [](auto& a, auto& b){ return a.first < b.first; });

    vector<pair<string,double>> elems(Imap.begin(), Imap.end());
    sort(elems.begin(), elems.end(), [](auto& a, auto& b){ return a.first < b.first; });

    bool quit=false; SDL_Event e; int scrollL=0, scrollR=0; const int pad = 16;

    while(!quit){
        while(SDL_PollEvent(&e)){
            if(e.type==SDL_QUIT) quit=true;
            else if(e.type==SDL_KEYDOWN && e.key.keysym.sym==SDLK_ESCAPE) quit=true;
            else if(e.type==SDL_MOUSEWHEEL){
                // شیفت = اسکرول ستون راست
                if ((SDL_GetModState() & KMOD_SHIFT) != 0) scrollR += (e.wheel.y>0? 20 : -20);
                else                                       scrollL += (e.wheel.y>0? 20 : -20);
            }
        }

        int W,H; SDL_GetWindowSize(w,&W,&H);
        SDL_SetRenderDrawColor(r,255,255,255,255); SDL_RenderClear(r);

        // عنوان
        drawText("DC Operating Point", pad, 8);
        drawText("(Ground = " + gndName + ")", pad+180, 8);

        // جداکنندهٔ عمودی وسط
        int midX = W/2;
        SDL_SetRenderDrawColor(r, 0,0,0,255);
        SDL_RenderDrawLine(r, midX, 40, midX, H-20);

        // ستون چپ: Node Voltages
        drawText("Node", pad, 40);
        drawText("Voltage (V)", pad+140, 40);
        SDL_RenderDrawLine(r, pad, 60, midX - pad, 60);

        int yL = 70 + scrollL;
        for (auto& kv : nodes) {
            if (kv.first == gndName) continue; // 0V
            drawText(kv.first, pad, yL);
            drawText(formatWithSI(kv.second, "V"), pad + 140, yL);
            yL += 22;
        }

// ستون راست: Element Currents
        drawText("Element", midX + pad, 40);
        drawText("Current (A)", midX + pad + 160, 40);
        SDL_RenderDrawLine(r, midX + pad, 60, W - pad, 60);

        int yR = 70 + scrollR;
        for (auto& kv : elems) {
            const string& ename = kv.first;
            double Iraw = kv.second;

            bool isNaN = isnan(Iraw);
            drawText(ename, midX + pad, yR);
            drawText(isNaN ? "N/A" : formatWithSI(Iraw, "A"), midX + pad + 160, yR);
            yR += 22;
        }

        SDL_RenderPresent(r);
    }

    if (f) TTF_CloseFont(f);
    SDL_DestroyRenderer(r);
    SDL_DestroyWindow(w);
}


void runCircuit() {
    if (!isGroundPlaced) return;

    buildCircuit(circuit);

    // --- جدید: محاسبه و نمایش DC‑OP (ولتاژها + جریان‌ها در یک پنجره) ---
    try {
        map<string,double> Vmap, Imap;
        circuit.computeDCOP(Vmap, Imap);

        // چاپ ترمینال (اختیاری)
        cout << "\nNodal Voltages (ground " << circuit.groundName << " = 0 V):\n";
        for (auto& kv : Vmap) cout << kv.first << " = " << kv.second << " V\n";
        cout << "\nElement Currents (A) [positive n1->n2]:\n";
        for (auto& kv : Imap) {
            if (isnan(kv.second)) cout << kv.first << " = N/A\n";
            else                        cout << kv.first << " = " << kv.second << " A\n";
        }

        showDCOPWindow(Vmap, Imap, circuit.groundName);
    } catch (const exception& ex) {
        cout << ex.what() << endl;
        return;
    }

    cout << "[Run]" << endl;
}

static vector<double> logspace(double fStart, double fStop, int pointsPerDec) {
    if (fStart <= 0.0 || fStop <= 0.0 || fStop < fStart) {
        throw runtime_error("Invalid AC sweep range.");
    }
    if (pointsPerDec < 1) pointsPerDec = 1;

    const double decades = log10(fStop) - log10(fStart);
    const int N = max(2, int(ceil(decades * pointsPerDec)) + 1);

    vector<double> f; f.reserve(N);
    const double logStart = log10(fStart);
    const double step     = (log10(fStop) - logStart) / double(N - 1);
    for (int i = 0; i < N; ++i) f.push_back(pow(10.0, logStart + i * step));
    // تضمین اولین/آخرین
    f.front() = fStart; f.back() = fStop;
    return f;
}

// حل AC در یک فرکانس با MNA (R/C/L + Vsrc/Isrc مستقل)
// نکته‌ها:
//  - در AC، DC-offset منابع نادیده گرفته می‌شود؛ از amplitude برای Vsrc استفاده شده.
//  - برای Isrc اگر AC-amplitude ندارید، فعلاً

static inline bool isValidProbeChar(char c) {
    return isalnum((unsigned char)c) || c == '_' || c == '-';
}

static string chooseProbeNode(const Circuit& C) {
    // اول: پایانه‌ی غیرزمینِ اولین منبع ولتاژ
    for (auto* e : C.elements) {  // ← توجه: بدون ()
        if (auto* v = dynamic_cast<VoltageSource*>(e)) {
            const string n1 = v->getNode1()->getName();
            const string n2 = v->getNode2()->getName();
            if (n1 != C.groundName) return n1;
            if (n2 != C.groundName) return n2;
        }
    }
    // اگر نبود: اولین نود غیرزمین
    for (auto* n : C.nodes) {      // ← توجه: بدون ()
        if (n->getName() != C.groundName) return n->getName();
    }
    return {};
}



static bool solveACAtFrequency(
        const Circuit& C, double freq,
        map<string, cd>& outV)
{
    outV.clear();
    if (freq <= 0.0) return false;


    // 1) فهرست گره‌های مجهول (به‌جز زمین)
    vector<Node*> nodes;
    nodes.reserve(C.nodes.size());
    for (auto* n : C.nodes) if (n->getName() != C.groundName) nodes.push_back(n);
    const int N = (int)nodes.size();
    if (N == 0) return true;

    // 2) ولتاژسورس‌ها
    vector<const VoltageSource*> vsrcs;
    vsrcs.reserve(C.elements.size());
    for (auto* e : C.elements) if (auto* v = dynamic_cast<VoltageSource*>(e)) vsrcs.push_back(v);
    const int M = (int)vsrcs.size();

    // 3) اندیس گره‌ها
    map<string,int> idx;
    for (int i=0;i<N;++i) idx[nodes[i]->getName()] = i;

    // 4) ماتریس MNA: [Y  B; B^T 0] x [V; Ivs] = [I; E]
    const int DIM = N + M;
    const double w = 2.0 * M_PI * freq;

    vector<vector<cd>> A(DIM, vector<cd>(DIM, cd(0,0)));
    vector<cd> rhs(DIM, cd(0,0));

    auto add = [&](int r, int c, cd v){ if (r>=0 && c>=0) A[r][c] += v; };

    // 5) استمپ R/C/L و Isrc
    for (auto* e : C.elements) {
        const string a = e->getNode1()->getName();
        const string b = e->getNode2()->getName();
        const int ia = (a != C.groundName) ? idx[a] : -1;
        const int ib = (b != C.groundName) ? idx[b] : -1;

        if (auto* R = dynamic_cast<Resistor*>(e)) {
            double r = R->getValue(); r = (r <= 1e-12 ? 1e-12 : r);
            cd y = cd(1.0 / r, 0.0);
            if (ia!=-1) add(ia, ia,  y);
            if (ib!=-1) add(ib, ib,  y);
            if (ia!=-1 && ib!=-1) { add(ia, ib, -y); add(ib, ia, -y); }
        }
        else if (auto* Cc = dynamic_cast<Capacitor*>(e)) {
            double Cv = Cc->getValue(); Cv = (Cv < 0 ? 0 : Cv);
            cd y = j * (w * Cv);
            if (ia!=-1) add(ia, ia,  y);
            if (ib!=-1) add(ib, ib,  y);
            if (ia!=-1 && ib!=-1) { add(ia, ib, -y); add(ib, ia, -y); }
        }
        else if (auto* L = dynamic_cast<Inductor*>(e)) {
            double Lv = L->getValue(); Lv = (Lv <= 0 ? 1e-12 : Lv);
            cd y = cd(1.0, 0.0) / (j * (w * Lv));
            if (ia!=-1) add(ia, ia,  y);
            if (ib!=-1) add(ib, ib,  y);
            if (ia!=-1 && ib!=-1) { add(ia, ib, -y); add(ib, ia, -y); }
        }
        else if (auto* I = dynamic_cast<CurrentSource*>(e)) {
            // جریان مثبت از n1→n2
            // رفتار اسپایسی: AC=0 مگر amplitude جداگانه داشته باشد.
            cd Iac = cd(0.0, 0.0);
            // اگر فیلدی مثل I->amplitude داری، از آن استفاده کن:
            // Iac = cd(I->amplitude, 0.0);
            if (ia!=-1) rhs[ia] +=  Iac;
            if (ib!=-1) rhs[ib] += -Iac;
        }
    }

    // 6) ولتاژسورس‌ها: ستون‌های اضافه و RHS
    for (int k=0; k<M; ++k) {
        const auto* vs = vsrcs[k];
        const string a = vs->getNode1()->getName();
        const string b = vs->getNode2()->getName();
        const int ia = (a != C.groundName) ? idx[a] : -1;
        const int ib = (b != C.groundName) ? idx[b] : -1;

        const int col = N + k;
        if (ia!=-1) { add(ia, col, cd(1,0)); add(col, ia, cd(1,0)); }
        if (ib!=-1) { add(ib, col, -cd(1,0)); add(col, ib, -cd(1,0)); }

        // RHS: در AC از amplitude استفاده کن؛ offset بی‌اثر است.
        rhs[col] += cd(vs->amplitude, 0.0);
        // اگر فاز خواستی بعداً: rhs[col] = amp * exp(j*phi);
    }

    // 7) حل دستگاه مختلط با گاوس (بدون تبدیل نادرست complex→double)
    auto solveComplex = [&](vector<vector<cd>> M, vector<cd> b)->vector<cd>{
        const int n = (int)M.size();
        for (int i=0;i<n;++i){
            // انتخاب پیوت
            int piv = i;
            for (int r=i+1;r<n;++r)
                if (abs(M[r][i]) > abs(M[piv][i])) piv = r;
            swap(M[i], M[piv]);
            swap(b[i], b[piv]);

            cd p = M[i][i];
            if (abs(p) < 1e-18) throw SingularMatrixException();
            for (int c=i;c<n;++c) M[i][c] /= p;
            b[i] /= p;
            for (int r=i+1;r<n;++r){
                cd f = M[r][i];
                if (f == cd(0,0)) continue;
                for (int c=i;c<n;++c) M[r][c] -= f * M[i][c];
                b[r] -= f * b[i];
            }
        }
        vector<cd> x(n, cd(0,0));
        for (int i=n-1;i>=0;--i){
            x[i] = b[i];
            for (int j=i+1;j<n;++j) x[i] -= M[i][j]*x[j];
        }
        return x;
    };


    (void)DIM; // سکوت هشدار unused

    vector<cd> X = solveComplex(A, rhs);

    // 8) خروجی: ولتاژ گره‌ها
    for (int i=0;i<N;++i) outV[nodes[i]->getName()] = X[i];
    return true;
}



// نمایش نتایج AC به‌صورت جدول (Freq, |V|, ∠V) شبیه پنجره DC
static void showACSweepWindow(const vector<double>& f,
                              const vector<double>& mag,
                              const vector<double>& pha_deg,
                              const string& nodeName)
{
    SDL_Window* w = SDL_CreateWindow(("AC Sweep: V(" + nodeName + ")").c_str(),
                                     SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                     720, 520, SDL_WINDOW_SHOWN);
    SDL_Renderer* r = SDL_CreateRenderer(w, -1, SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* fnt = TTF_OpenFont("C:\\Windows\\Fonts\\Arial.ttf", 14);

    auto text = [&](const string& s, int x, int y){
        if (!fnt) return;
        SDL_Color c{0,0,0};
        SDL_Surface* sf = TTF_RenderText_Solid(fnt, s.c_str(), c);
        if (!sf) return;
        SDL_Texture* tx = SDL_CreateTextureFromSurface(r, sf);
        SDL_Rect dst{x, y, sf->w, sf->h};
        SDL_RenderCopy(r, tx, nullptr, &dst);
        SDL_DestroyTexture(tx); SDL_FreeSurface(sf);
    };

    int scroll=0; bool quit=false; SDL_Event e;
    while(!quit){
        while(SDL_PollEvent(&e)){
            if (e.type==SDL_QUIT) quit=true;
            else if (e.type==SDL_KEYDOWN && e.key.keysym.sym==SDLK_ESCAPE) quit=true;
            else if (e.type==SDL_MOUSEWHEEL) scroll += (e.wheel.y>0? 20 : -20);
        }
        SDL_SetRenderDrawColor(r,255,255,255,255); SDL_RenderClear(r);
        text("Freq(Hz)", 20, 10);  text("|V| (V)", 200, 10);  text("Phase (deg)", 330, 10);
        SDL_SetRenderDrawColor(r,0,0,0,255);
        SDL_RenderDrawLine(r, 10, 30, 700, 30);

        int y = 50 + scroll;
        for (size_t i=0;i<f.size();++i){
            text(fmt6(f[i]),          20,  y);
            text(fmt6(mag[i]),       200,  y);
            text(fmt6(pha_deg[i]),   330,  y);
            y += 22;
        }
        SDL_RenderPresent(r);
    }
    if (fnt) TTF_CloseFont(fnt);
    SDL_DestroyRenderer(r); SDL_DestroyWindow(w);
}


static vector<double> logspace(double fStart, double fStop, int num); // همونی که داری/یا معادل

// اگر probeName خالی باشد، مثل قبل پروب خودکار انتخاب می‌شود
void runACSweepFromCanvas(Circuit& C,
                          double fStart, double fStop, int ptsPerDec,
                          const string& probeName /* may be "" */)
{
    buildCircuit(C);

    if (!(fStart > 0 && fStop > 0 && fStop >= fStart)) {
        cout << "AC sweep params invalid.\n"; return;
    }
    if (ptsPerDec < 1) ptsPerDec = 1;

    if (C.groundName.empty()) {
        cout << "Ground node is not set.\n"; return;
    }

    // انتخاب پروب: یا از ورودی کاربر یا خودکار
    string probe = probeName;
    if (!probe.empty()) {
        // ولیدیشن: گره وجود دارد و زمین نیست
        bool exists = false;
        for (auto* n : C.nodes) if (n->getName() == probe) { exists = true; break; }
        if (!exists) { cout << "Probe node not found: " << probe << "\n"; return; }
        if (probe == C.groundName) { cout << "Probe cannot be ground.\n"; return; }
    } else {
        probe = chooseProbeNode(C);
        if (probe.empty()) { cout << "No non-ground node to probe.\n"; return; }
    }

    // تولید فرکانس‌ها (logspace مشابه قبل)
    const double decades = log10(fStop / fStart);
    const int num = max(2, int(round(decades * ptsPerDec)) + 1);
    auto freqs = logspace(fStart, fStop, num);

    vector<double> mags; mags.reserve(freqs.size());
    vector<double> phs;  phs.reserve(freqs.size());

    for (double f : freqs) {
        map<string, cd> V;
        if (!solveACAtFrequency(C, f, V)) { mags.push_back(0.0); phs.push_back(0.0); continue; }
        auto it = V.find(probe);
        cd v = (it == V.end()? cd(0.0,0.0) : it->second);
        mags.push_back(abs(v));
        phs.push_back(arg(v) * 180.0 / M_PI);
    }

    showACSweepWindow(freqs, mags, phs, "V(" + probe + ")");
}




// State for AC Sweep dialog
string acStartText = "10";    // Hz
string acStopText  = "10k";   // Hz
string acPointsPerDecText = "10"; // points/decade
string acProbeText = ""; // points/decade
int acFocus = -1; // 0=start, 1=stop, 2=points


// ===== Run dialog state =====
bool showRunDialog = false;

enum RunChoice { RUN_DC, RUN_VT, RUN_IT , RUN_AC };
RunChoice runChoice = RUN_DC;

string runNodeName = "";   // برای V–t
string runElemName = "";   // برای I–t

enum RunField { RF_NONE, RF_NODE, RF_ELEM };
RunField runFocus = RF_NONE;



SDL_Rect runDialogRect = { 360, 180, 560, 200 }; // جای دیالوگ Run

// اعلان توابع نمایش سیگنال
void showSignalVT(const vector<double>& t,
                  const vector<double>& y,
                  const string& title,
                  const string& yLabel);
void showNodeVoltageVT(const WaveStore& ws, const string& nodeName);
void showElementCurrentIT(const WaveStore& ws, const string& elemName);


// اگر هنوز simulateTransientCapture نداری، تابعش در بخش 2 آمده.
void showSignalVT(const vector<double>& t,
                  const vector<double>& y,
                  const string& title,
                  const string& yLabel) {
    if (t.empty() || y.empty() || t.size() != y.size()) return;

    SDL_Window* w = SDL_CreateWindow(title.c_str(),
                                     SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 900, 480, SDL_WINDOW_SHOWN);
    SDL_Renderer* r = SDL_CreateRenderer(w, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    TTF_Font* font2 = TTF_OpenFont("C:\\Windows\\Fonts\\Arial.ttf", 14);

    double tMin = t.front(), tMax = t.back();
    double yMin = y[0], yMax = y[0];
    for (size_t i = 1; i < y.size(); ++i) { yMin = min(yMin, y[i]); yMax = max(yMax, y[i]); }
    if (fabs(yMax - yMin) < 1e-12) { yMax = yMin + 1.0; }

    auto X = [&](double tt, int W, int pad){ return pad + int((tt - tMin) * (W - 2*pad) / (tMax - tMin)); };
    auto Y = [&](double yy, int H, int pad){ return H - pad - int((yy - yMin) * (H - 2*pad) / (yMax - yMin)); };

    bool quit = false;
    SDL_Event e;
    const int pad = 50;

    auto drawText = [&](const string& s, int x, int yPos){
        if (!font2) return;
        SDL_Color black = {0,0,0};
        SDL_Surface* surf = TTF_RenderText_Solid(font2, s.c_str(), black);
        if (!surf) return;
        SDL_Texture* tex = SDL_CreateTextureFromSurface(r, surf);
        SDL_Rect dst = {x, yPos, surf->w, surf->h};
        SDL_RenderCopy(r, tex, nullptr, &dst);
        SDL_DestroyTexture(tex); SDL_FreeSurface(surf);
    };

    while (!quit) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = true;
            else if (e.type == SDL_KEYDOWN && (e.key.keysym.sym == SDLK_ESCAPE)) quit = true;
        }

        int W, H; SDL_GetWindowSize(w, &W, &H);
        SDL_SetRenderDrawColor(r, 255, 255, 255, 255); SDL_RenderClear(r);

        // axes
        SDL_SetRenderDrawColor(r, 0, 0, 0, 255);
        SDL_RenderDrawLine(r, pad, H - pad, W - pad, H - pad); // t-axis
        SDL_RenderDrawLine(r, pad, pad,     pad,     H - pad); // y-axis

        drawText("t (s)", W - pad - 40, H - pad + 8);
        drawText(yLabel, 10, 10);
        drawText(title, W/2 - 80, 10);

        // line
        for (size_t i = 1; i < t.size(); ++i) {
            int x1 = X(t[i-1], W, pad), y1 = Y(y[i-1], H, pad);
            int x2 = X(t[i],   W, pad), y2 = Y(y[i],   H, pad);
            SDL_RenderDrawLine(r, x1, y1, x2, y2);
        }

        SDL_RenderPresent(r);
    }

    if (font2) TTF_CloseFont(font2);
    SDL_DestroyRenderer(r);
    SDL_DestroyWindow(w);
}

void showNodeVoltageVT(const WaveStore& ws, const string& nodeName) {
    auto it = ws.V.find(nodeName);
    if (it == ws.V.end()) return;
    showSignalVT(ws.t, it->second, "V(" + nodeName + ") vs t", "V (V)");
}

void showElementCurrentIT(const WaveStore& ws, const string& elemName) {
    auto it = ws.I.find(elemName);
    if (it == ws.I.end()) return;
    showSignalVT(ws.t, it->second, "I(" + elemName + ") vs t", "I (A)");
}




// گرد کردن به نزدیک‌ترین گره شبکه
inline void snapToGrid(int x, int y, int& gx, int& gy) {
    gx = (int)round((double)x / GRID_SPACING) * GRID_SPACING;
    gy = (int)round((double)y / GRID_SPACING) * GRID_SPACING;
}

// اولویت: کانکتورهای موجود → شبکه
inline void snapPointWithConnectorsFirst(int x, int y, int& outX, int& outY) {
    // 1) اگر نزدیک کانکتور موجود باشیم، به همان اسنپ کن
    for (const auto& c : allConnectors) {
        if (abs(x - c.x) <= SNAP_BOX_HALF && abs(y - c.y) <= SNAP_BOX_HALF) {
            outX = c.x; outY = c.y;
            return;
        }
    }
    // 2) در غیر اینصورت، به شبکه اسنپ کن (در محدوده جذب 10x10)
    int gx, gy; snapToGrid(x, y, gx, gy);
    if (abs(x - gx) <= SNAP_BOX_HALF && abs(y - gy) <= SNAP_BOX_HALF) {
        outX = gx; outY = gy;
    } else {
        // خارج محدوده جذب: اجازه نده روی بوم چیزی گذاشته شود (می‌تونی به نقطه خام برگردونی)
        outX = gx; outY = gy; // پیشنهاد: همچنان به نزدیک‌ترین گره بچسبانیم
    }
}

inline void drawGrid(SDL_Renderer* r, int W, int H) {
    SDL_SetRenderDrawColor(r, 220, 220, 220, 255); // خاکستری روشن
    for (int y = 0; y <= H; y += GRID_SPACING) {
        for (int x = 0; x <= W; x += GRID_SPACING) {
            SDL_Rect dot = { x - GRID_DOT_SIZE/2, y - GRID_DOT_SIZE/2, GRID_DOT_SIZE, GRID_DOT_SIZE };
            SDL_RenderFillRect(r, &dot);
        }
    }
}



// main headers

SCircuit scircuit;

string path = "";

struct XYKey { int x,y; };
struct XYHash {
    size_t operator()(const XYKey& k) const noexcept {
        return (uint64_t(uint32_t(k.x))<<32) ^ uint32_t(k.y);
    }
};
struct XYEq {
    bool operator()(const XYKey& a, const XYKey& b) const noexcept {
        return a.x==b.x && a.y==b.y;
    }
};

// گرفتن/ساختن connector id برای یک نقطه
inline int get_or_make_connector(int x, int y,
                                 std::unordered_map<XYKey,int,XYHash,XYEq>& xy2cid)
{
    XYKey key{x,y};
    auto it = xy2cid.find(key);
    if (it != xy2cid.end()) return it->second;

    const int newId = (int)allConnectors.size() + 1; // یا از connectorIdCounter++
    allConnectors.push_back(Connector{ newId, x, y });
    xy2cid.emplace(key, newId);
    return newId;
}


//enum ToolType { RESISTOR, CAPACITOR, INDUCTOR, VSOURCE, CSOURCE, GROUND, VSIN, WIRE };

static ToolType kind_to_tool(const std::string& k){
    if (k=="R")   return RESISTOR;
    if (k=="C")   return CAPACITOR;
    if (k=="L")   return INDUCTOR;
    if (k=="V")   return VSOURCE; // VSIN را هم در value نگه داشتیم؛ اگر خواستی جداگانه VSIN بده
    if (k=="I")   return CSOURCE;
    if (k=="GND") return GROUND;
    return RESISTOR;
}


// این تابع: فایل را باز می‌کند، SCircuit می‌خواند، و بوم را می‌سازد
bool openAndLoadToCanvas()
{
    // 1) مسیر فایل
    auto path = show_open_dialog_win("Open Circuit (.json)", "json");
    if (path.empty()) return false;

    // 2) خواندن SCircuit
    SCircuit ckt;
    if (!loadProject(path, ckt)) return false;

    // 3) پاک کردن بوم فعلی
    placedElements.clear();
    wires.clear();
    allConnectors.clear();
    connectorToNode.clear();

    // 4) map برای تبدیل node_id → name  (و برعکس اگر خواستی)
    std::vector<std::string> nodeIdToName;
    nodeIdToName.resize(ckt.Snodes.size());
    for (const auto& n : ckt.Snodes) {
        if (n.id >= 0 && (size_t)n.id < nodeIdToName.size())
            nodeIdToName[n.id] = n.name;
    }

    // 5) نگاشت مختصات → connector id
    std::unordered_map<XYKey,int,XYHash,XYEq> xy2cid; xy2cid.reserve(256);

    // 6) بازسازی المان‌ها و کانکتورهایشان
    const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB;

    for (const auto& se : ckt.Selements) {
        // تبدیل kind → ToolType (ساده)
        ToolType tt = kind_to_tool(se.kind);

        // ساخت PlacedElement (بدون چرخش/وارونه)
        PlacedElement pe{};
        pe.type = tt;
        pe.x    = se.x;
        pe.y    = se.y;
        pe.value= se.value;
        pe.amp  = 0;          // اگر برای VSIN بعداً نیاز شد پر کن
        pe.freq = 0;
        pe.name = se.name;

        // اضافه به بوم
        placedElements.push_back(pe);

        // محاسبهٔ کانکتورهای المان (افقی ثابت)
        int id1, id2;
        if (tt == GROUND) {
            id1 = get_or_make_connector(se.x, se.y, xy2cid);
            id2 = id1;
        } else {
            id1 = get_or_make_connector(se.x - PIN_OFFSET, se.y, xy2cid);
            id2 = get_or_make_connector(se.x + PIN_OFFSET, se.y, xy2cid);
        }

        // نگاشت connector → node name  (برای هر پین)
        if (se.n1 >= 0 && (size_t)se.n1 < nodeIdToName.size()) {
            connectorToNode[id1] = nodeIdToName[se.n1];
            loadedSeed[id1]      = nodeIdToName[se.n1];
        }
        if (se.n2 >= 0 && (size_t)se.n2 < nodeIdToName.size()) {
            connectorToNode[id2] = nodeIdToName[se.n2];
            loadedSeed[id2]      = nodeIdToName[se.n2];
        }


        if (tt == GROUND) {

            int gnid = (se.n1 >= 0) ? se.n1 : se.n2;
            if (gnid >= 0 && (size_t)gnid < nodeIdToName.size()) {
                const std::string& gname = nodeIdToName[gnid];

                // کانکتور نماد زمین را حتماً به همون نود مپ کن
                connectorToNode[id1] = gname;

                // فلگ UI تا اجرای Run بدون خطا
                isGroundPlaced = true;

                if (circuit.getGroundName().empty()) circuit.setGroundNode(gname);
            }
        }
    }

    // 7) بازسازی سیم‌ها (Polyline → سگمنت‌های متوالی)
    int nextWireId = 1;
    for (const auto& sw : ckt.Swires) {
        if (sw.points.size() < 2) continue;
        for (size_t i = 1; i < sw.points.size(); ++i) {
            const auto& a = sw.points[i-1];
            const auto& b = sw.points[i];
            // حذف صفرطول
            if (a.x==b.x && a.y==b.y) continue;

            // تبدیل به WireSeg برنامه‌ی تو
            SDL_Point p1{ a.x, a.y };
            SDL_Point p2{ b.x, b.y };
            wires.push_back({ p1, p2 });


            // کانکتورهای متناظر هر سر
            int c1 = get_or_make_connector(a.x, a.y, xy2cid);
            int c2 = get_or_make_connector(b.x, b.y, xy2cid);

            // اگر برای سر سیم‌ها node داریم، می‌تونی اینجا connectorToNode رو هم پر کنی؛
            // اما چون node سیم از روی اتصال با المان‌ها/برچسب‌ها مشخص میشه، معمولا همین کافیست.
            (void)c1; (void)c2;
            ++nextWireId;
        }
    }

    // 8) GND به‌عنوان عنصر: اگر می‌خواهی همان‌جا ground را ست کنی
    // هنگام تحلیل بعدی، اگر kind=="GND"، آن نود را زمین در نظر بگیر.

    // 9) موفق
    return true;
}








bool saveFromCanvasLikeBuildCircuit(const std::string& folderPath)
{
    // 0) بازسازی توپولوژی (عین buildCircuit)
    allConnectors.clear();
    connectorToNode.clear();
    analyzeNodeConnections();

    // 1) مختصات نماینده برای هر نود منطقی
    std::unordered_map<std::string, std::pair<int,int>> nodePos;
    nodePos.reserve(allConnectors.size());
    for (const auto& c : allConnectors) {
        auto it = connectorToNode.find(c.id);
        if (it == connectorToNode.end()) continue;
        const std::string& nn = it->second;
        if (!nodePos.count(nn)) nodePos[nn] = {c.x, c.y};
    }

    // 2) ساخت SCircuit و Snodes
    SCircuit c;
    c.Snodes.reserve(nodePos.size());
    std::unordered_map<std::string,int> nodeId; nodeId.reserve(nodePos.size());
    int nextId = 0;
    for (auto& kv : nodePos) {
        const std::string& name = kv.first;
        auto xy = kv.second;
        nodeId[name] = nextId;
        c.Snodes.push_back(SNode{ nextId, name, xy.first, xy.second });
        ++nextId;
    }

    // 3) نام‌گذاری خودکار مثل buildCircuit
    const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB;
    std::map<ToolType,int> typeCounter;
    auto autoName = [&](const PlacedElement& e)->std::string{
        static std::map<ToolType,char> L = {
                {RESISTOR,'R'},{CAPACITOR,'C'},{INDUCTOR,'L'},
                {VSOURCE,'V'},{CSOURCE,'I'},{GROUND,'G'},{VSIN,'V'}
        };
        if (!e.name.empty()) return e.name;
        return std::string(1, L[e.type]) + std::to_string(++typeCounter[e.type]);
    };

    // 4) تبدیل placedElements → Selements (بدون rot/mir)
    c.Selements.reserve(placedElements.size());
    for (const auto& e : placedElements) {
        if (e.type == WIRE) continue;

        int id1, id2;
        const int y = e.y;
        if (e.type == GROUND) {
            id1 = addConnector(e.x, y);
            id2 = id1;
        } else {
            id1 = addConnector(e.x - PIN_OFFSET, y);
            id2 = addConnector(e.x + PIN_OFFSET, y);
        }

        auto it1 = connectorToNode.find(id1);
        auto it2 = connectorToNode.find(id2);
        if (it1 == connectorToNode.end() || it2 == connectorToNode.end()) continue;

        const std::string node1 = it1->second;
        const std::string node2 = it2->second;

        int n1 = -1, n2 = -1;
        auto f1 = nodeId.find(node1); if (f1 != nodeId.end()) n1 = f1->second;
        auto f2 = nodeId.find(node2); if (f2 != nodeId.end()) n2 = f2->second;

        std::string name = autoName(e);

        std::string kind;
        switch (e.type) {
            case RESISTOR: kind = "R"; break;
            case CAPACITOR: kind = "C"; break;
            case INDUCTOR: kind = "L"; break;
            case VSOURCE:  kind = "V"; break;   // DC
            case VSIN:     kind = "V"; break;   // AC (فعلاً فقط value=offset)
            case CSOURCE:  kind = "I"; break;
            case GROUND:   kind = "GND"; break;
            default: continue;
        }

        // چرخش/آینه نداریم: rotation=0.0, mirror=false
        c.Selements.push_back(SElement{
                kind, name,
                n1, n2,
                e.value,
                e.x, e.y,
                0.0, false
        });
    }

    // 5) ذخیره سیم‌ها به‌صورت گرافیکی (Polyline دو نقطه‌ای)
    c.Swires.reserve(wires.size());
    for (const auto& w : wires) {
        SWire sw; sw.id = static_cast<int>(c.Swires.size()) + 1;
        sw.points.push_back(Point{ w.first.x,  w.first.y  });
        sw.points.push_back(Point{ w.second.x, w.second.y });
        c.Swires.push_back(std::move(sw));
    }

    // 6) ذخیره روی دیسک
    return saveProject(folderPath, c);
}


//================ Main Function =============================
int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    TTF_Init();
    SDL_Window* window = SDL_CreateWindow("Opulator Spice", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    TTF_Font* font = TTF_OpenFont("C:\\Windows\\Fonts\\Arial.ttf", 14);

    int winW, winH; SDL_GetWindowSize(window, &winW, &winH);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
//    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    SDL_RenderPresent(renderer);

    bool showFileMenu = false;
    bool fileClicked = false;

    bool showPreview = false;
    int preX, preY;

    bool showValueDialog = false;
    ToolType selectedToolForDialog = NONE;
    string inputValueText = "";
    SDL_Rect valueDialogRect = { 400, 300, 300, 120 };
    SDL_Rect vacDialogRect = { 380, 220, 420, 260 };
    double pendingValue = 0.0;   // مقدار بعد از OK
    bool   hasPendingValue = false;
    bool showVACDialog = false;
    string vacNameText, vacAmpText, vacFreqText, vacPhaseText;
    int vacFocus = 0;
    struct PendingVAC { string name; double  offset = 0.0,amp=0, freq=0, phase=0; } pendingVAC;
    string vacOffsetText;
    // ===== Dialog state (add these) =====
    enum DialogField { FIELD_NAME, FIELD_VALUE };
    DialogField dialogFocus = FIELD_NAME;  // فوکوس اول روی نام

    string inputNameText = "";        // ورودی نام

    string pendingName;               // نام تایید شده بعد از OK
    bool hasPendingName = false;    // فلگ در دسترس بودن نام

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
            }
            else if (e.type == SDL_MOUSEMOTION) {
                showFileMenu = hoveringFile || (hoveringDropdown && fileClicked);

                if (currentTool != NONE && e.motion.y > 75) {
                    int sx, sy;
                    snapPointWithConnectorsFirst(e.motion.x, e.motion.y, sx, sy);
                    preX = sx; preY = sy;
                    showPreview = true;

                }

                else {
                    showPreview = false;
                }
            }


            else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                cout << "x: " << e.button.x << " y: " << e.button.y << endl;
                // --- RUN DIALOG CLICK HANDLER ---
                // --- FIX 1: put this lambda once near start of event processing (inside while(SDL_PollEvent))
                auto hit = [&](SDL_Rect r, int mx, int my)->bool {
                    return mx >= r.x && mx <= r.x + r.w && my >= r.y && my <= r.y + r.h;
                };

// ... inside the existing mouse-left handler (no second vacDialogRect definition!)

                if (showVACDialog) {
                    SDL_Rect okBtn     = { vacDialogRect.x + vacDialogRect.w - 210, vacDialogRect.y + vacDialogRect.h - 50, 90, 28 };
                    SDL_Rect cancelBtn = { okBtn.x + 105, okBtn.y, 90, 28 };
                    int mx = e.button.x, my = e.button.y;

                    if (hit(okBtn, mx, my)) {
                        try {
                            pendingVAC.name   = vacNameText.empty()? "VAC" : vacNameText;
                            pendingVAC.offset = vacOffsetText.empty()? 0.0 : parseNumber(vacOffsetText);
                            pendingVAC.amp    = parseNumber(vacAmpText);
                            pendingVAC.freq   = parseNumber(vacFreqText);
                            pendingVAC.phase  = vacPhaseText.empty()? 0.0 : parseNumber(vacPhaseText);

                            showVACDialog = false;
                            SDL_StopTextInput();
                            currentTool = VSIN;   // وارد پیش‌نمایش
                            showPreview = true;
                            continue; // جلوگیری از نفوذ
                        } catch(...) {
                            cout << "Invalid AC inputs!\n";
                        }
                    } else if (hit(cancelBtn, mx, my)) {
                        showVACDialog = false;
                        SDL_StopTextInput();
                        continue;
                    }

                    // فوکوس با کلیک روی باکس‌ها (اختیاری):
                    int xL = vacDialogRect.x + 14, xR = vacDialogRect.x + 220;
                    int y0 = vacDialogRect.y + 16, gap = 42, boxW = 170, boxH = 28;
                    SDL_Rect boxName  { xR, y0-4,        boxW, boxH };
                    SDL_Rect boxOffs  { xR, y0+gap-4,    boxW, boxH };
                    SDL_Rect boxAmp   { xR, y0+2*gap-4,  boxW, boxH };
                    SDL_Rect boxFreq  { xR, y0+3*gap-4,  boxW, boxH };
                    SDL_Rect boxPhase { xR, y0+4*gap-4,  boxW, boxH };
                    if      (hit(boxName,  mx,my)) { vacFocus = 0; SDL_StartTextInput(); continue; }
                    else if (hit(boxOffs,  mx,my)) { vacFocus = 1; SDL_StartTextInput(); continue; }
                    else if (hit(boxAmp,   mx,my)) { vacFocus = 2; SDL_StartTextInput(); continue; }
                    else if (hit(boxFreq,  mx,my)) { vacFocus = 3; SDL_StartTextInput(); continue; }
                    else if (hit(boxPhase, mx,my)) { vacFocus = 4; SDL_StartTextInput(); continue; }

                    continue; // وقتی دیالوگ باز است، کلیک به پایین دسترس نداشته باشد
                }

                if (showVACDialog) {
                    SDL_Rect okBtn     = { vacDialogRect.x + vacDialogRect.w - 210, vacDialogRect.y + vacDialogRect.h - 40, 90, 28 };
                    SDL_Rect cancelBtn = { vacDialogRect.x + vacDialogRect.w - 105, vacDialogRect.y + vacDialogRect.h - 40, 90, 28 };
                    int mx = e.button.x, my = e.button.y;

                    if (mx >= okBtn.x && mx <= okBtn.x + okBtn.w && my >= okBtn.y && my <= okBtn.y + okBtn.h) {
                        try {
                            pendingVAC.name  = vacNameText.empty()? "VAC" : vacNameText;
                            pendingVAC.amp   = parseNumber(vacAmpText);    // مثل فاز۱
                            pendingVAC.freq  = parseNumber(vacFreqText);
                            pendingVAC.phase = vacPhaseText.empty()? 0.0 : parseNumber(vacPhaseText);

                            showVACDialog = false;
                            SDL_StopTextInput();
                            currentTool = VSIN;   // وارد فاز پیش‌نمایش
                            showPreview = true;
                            continue; // از نفوذ رویداد جلوگیری کن
                        } catch (...) {
                            cout << "Invalid AC inputs!" << endl;
                        }
                    } else if (mx >= cancelBtn.x && mx <= cancelBtn.x + cancelBtn.w &&
                               my >= cancelBtn.y && my <= cancelBtn.y + cancelBtn.h) {
                        showVACDialog = false;
                        SDL_StopTextInput();
                        continue;
                    }
                }
                if (showRunDialog) {
                    // مستطیل‌ها
                    SDL_Rect okBtn     = {runDialogRect.x + runDialogRect.w - 200, runDialogRect.y + runDialogRect.h - 40, 80, 28};
                    SDL_Rect cancelBtn = {okBtn.x + 100, okBtn.y, 80, 28};
                    SDL_Rect rbDC      = {runDialogRect.x + 20,  runDialogRect.y + 40,  18, 18};
                    SDL_Rect rbAC      = {runDialogRect.x + 20,  runDialogRect.y + 70,  18, 18};

                    // جعبه‌های AC
                    SDL_Rect acStartBox = {runDialogRect.x + 180, runDialogRect.y + 65,  80, 24};
                    SDL_Rect acStopBox  = {runDialogRect.x + 280, runDialogRect.y + 65,  80, 24};
                    SDL_Rect acPtsBox   = {runDialogRect.x + 380, runDialogRect.y + 65,  80,  24};
                    SDL_Rect acProbeBox = {runDialogRect.x + 480, runDialogRect.y + 65,  80,  24};

                    auto inRect = [&](SDL_Rect r)->bool{
                        return e.button.x >= r.x && e.button.x <= r.x + r.w &&
                               e.button.y >= r.y && e.button.y <= r.y + r.h;
                    };

                    // OK
                    if (inRect(okBtn)) {
                        SDL_StopTextInput();
                        if (runChoice == RUN_DC) {
                            showRunDialog = false;
                            runCircuit();
                        } else if (runChoice == RUN_AC) {
                            showRunDialog = false;
                            double fStart = parseNumber(acStartText);
                            double fStop  = parseNumber(acStopText);
                            int    pDec   = max(1, atoi(acPointsPerDecText.c_str()));
                            // قبلی:
                            // runACSweepStub(circuit, fStart, fStop, pDec);

                            // جدید:
                            runACSweepFromCanvas(circuit, fStart, fStop, pDec, acProbeText);

                        } else if (runChoice == RUN_VT) {
                            string node = runNodeName.empty()? "n1" : runNodeName;
                            showRunDialog = false;
                            circuit.simulateTransientCapture(gWaves);
                            showNodeVoltageVT(gWaves, node);
                        } else { // RUN_IT
                            string elem = runElemName.empty()? "R1" : runElemName;
                            showRunDialog = false;
                            circuit.simulateTransientCapture(gWaves);
                            showElementCurrentIT(gWaves, elem);
                        }
                        continue;
                    }

                    // Cancel
                    if (inRect(cancelBtn)) { showRunDialog = false; SDL_StopTextInput(); continue; }

                    // انتخاب نوع تحلیل
                    if (inRect(rbDC))  { runChoice = RUN_DC; acFocus = -1; runFocus = RF_NONE; continue; }
                    if (inRect(rbAC))  { runChoice = RUN_AC; runFocus = RF_NONE; acFocus = 0; SDL_StartTextInput(); continue; }

                    // فوکوس روی فیلدهای AC وقتی RUN_AC انتخاب شده
                    if (runChoice == RUN_AC) {
                        if (inRect(acStartBox)) { acFocus = 0; SDL_StartTextInput(); continue; }
                        if (inRect(acStopBox))  { acFocus = 1; SDL_StartTextInput(); continue; }
                        if (inRect(acPtsBox))   { acFocus = 2; SDL_StartTextInput(); continue; }
                        if (inRect(acProbeBox)) { acFocus = 3; SDL_StartTextInput(); continue; }
                    }

                    // کلیک بیرون از دیالوگ → نادیده
                    if (!(e.button.x >= runDialogRect.x && e.button.x <= runDialogRect.x + runDialogRect.w &&
                          e.button.y >= runDialogRect.y && e.button.y <= runDialogRect.y + runDialogRect.h)) {
                        continue;
                    }
                }
                SDL_Rect vacDialogRect = { 380, 220, 420, 220 }; // اگر بیرون اسکوپ لازمش نداری، همین‌جا تعریف باشه


                if (showValueDialog) {
                    // نواحی قابل کلیک داخل دیالوگ
                    SDL_Rect okButton = {valueDialogRect.x + 100, valueDialogRect.y + 80, 100, 30};
                    SDL_Rect nameBox  = {valueDialogRect.x + 60,  valueDialogRect.y + 8,  valueDialogRect.w - 70, 28};
                    SDL_Rect valBox   = {valueDialogRect.x + 60,  valueDialogRect.y + 38, valueDialogRect.w - 70, 28};

                    // اگر روی OK کلیک شد
                    if (e.button.x >= okButton.x && e.button.x <= okButton.x + okButton.w &&
                        e.button.y >= okButton.y && e.button.y <= okButton.y + okButton.h) {
                        try {
                            if (inputNameText.empty()) throw runtime_error("empty name");
                            double val = parseNumber(inputValueText);  // k/u/m/e… پشتیبانی می‌شود

                            preX = 545;
                            preY = 347;

                            pendingValue = val;
                            hasPendingValue = true;
                            pendingName = inputNameText;
                            hasPendingName = true;

                            showValueDialog = false;
                            SDL_StopTextInput();

                            currentTool = selectedToolForDialog;
                            selectedToolForDialog = NONE;
                            showPreview = true;

                            continue; // جلوگیری از نفوذ این کلیک به بقیه‌ی هندلرها
                        } catch (...) {
                            cout << "Invalid name/value!" << endl;
                        }
                        continue;
                    }

                    if (e.button.x >= nameBox.x && e.button.x <= nameBox.x + nameBox.w &&
                        e.button.y >= nameBox.y && e.button.y <= nameBox.y + nameBox.h) {
                        dialogFocus = FIELD_NAME;
                        continue;
                    }
                    if (e.button.x >= valBox.x && e.button.x <= valBox.x + valBox.w &&
                        e.button.y >= valBox.y && e.button.y <= valBox.y + valBox.h) {
                        dialogFocus = FIELD_VALUE;
                        continue;
                    }

                    continue;
                }


                if (hoveringDropdown && showFileMenu) {
                    int option = (my - 30) / 20;
                    switch (option) {
                        case 0:
                            placedElements.clear();
                            wires.clear();
                            allConnectors.clear();
                            connectorToNode.clear();

                            circuit.reset();
                            connectorIdCounter = 0;

                            path = "";
                            break;
                        case 1:
                            openAndLoadToCanvas();
                            break;
                        case 2:
                            if (path == "") {
                                path = show_save_dialog_win("Save Circuit", "circuit.json", "json");

                                if (path.empty()) {
                                    // کاربر لغو کرد
                                    continue;
                                }

                                saveFromCanvasLikeBuildCircuit(path);
                            } else {
                                saveFromCanvasLikeBuildCircuit(path);
                            }

                            break;
                    }
                    fileClicked = false;
                    showFileMenu = false;
                }  else if (mx >= 120 && mx <= 220 && my >= 0 && my <= 30) {
                    // اگر دیالوگ مقدار بازه، فعلاً کاری نکن تا تداخل پیش نیاد
                    if (!showValueDialog) {
                        showRunDialog = true;
                        runChoice = RUN_DC;
                        runNodeName.clear();
                        runElemName.clear();
                        runFocus = RF_NONE;
                        SDL_StartTextInput();
                    }
                }
                else if (my >= 25 && my <= 75) {
                    if (mx >= 30 && mx <= 60) currentTool = RESISTOR;
                    else if (mx >= 90 && mx <= 120) currentTool = CAPACITOR;
                    else if (mx >= 150 && mx <= 180) currentTool = INDUCTOR;
                    else if (mx >= 210 && mx <= 240) currentTool = VSOURCE;
                    else if (mx >= 270 && mx <= 300) currentTool = CSOURCE;
                    else if (mx >= 330 && mx <= 360) currentTool = GROUND;
                    else if (mx >= 390 && mx <= 420) currentTool = WIRE;
                    else if (mx >= 450 && mx <= 480) currentTool = VSIN;

                    // محدوده‌ی ابزارهای غیر از GROUND/WIRE: 30..330
                    if (currentTool == VSIN) {
                        showVACDialog = true;
                        vacNameText.clear();
                        vacOffsetText.clear();
                        vacAmpText.clear();
                        vacFreqText.clear();
                        vacPhaseText.clear();
                        vacFocus = 0;
                        currentTool = NONE;        // فعلاً ابزار غیرفعال تا OK
                        SDL_StartTextInput();
                    }

                    if (mx >= 30 && mx <= 330) {

                        if (currentTool != GROUND && currentTool != WIRE) {
                            selectedToolForDialog = currentTool;
                            currentTool = NONE;
                            showValueDialog = true;
                            inputNameText.clear();             // ← نام خالی
                            inputValueText.clear();            // ← مقدار خالی
                            dialogFocus = FIELD_NAME;          // ← فوکوس ابتدا روی نام
                            SDL_StartTextInput();              // ← شروع دریافت ورودی متن
                        }
                    }
                } else if (my > 75) {

                    if (currentTool == WIRE) {
                        if (!wireStartSelected) {
                            int sx, sy;
                            snapPointWithConnectorsFirst(mx, my, sx, sy);
                            wireX1 = sx; wireY1 = sy;
                            wireStartSelected = true;
                        } else {
                            int ex, ey;
                            snapPointWithConnectorsFirst(mx, my, ex, ey);
                            wires.push_back({ {wireX1, wireY1}, {ex, ey} });
                            wireStartSelected = false;
                            currentTool = NONE;
                            showPreview = false;
                        }
//                        if (!wireStartSelected) {
//                            // Snap start of wire to nearest connector
//                            int sx = mx, sy = my;
//                            for (auto& c : allConnectors) {
//                                if (abs(sx - c.x) <= 2 && abs(sy - c.y) <= 2) {
//                                    sx = c.x; sy = c.y;
//                                    break;
//                                }
//                            }
//                            wireX1 = sx; wireY1 = sy;
//                            wireStartSelected = true;
//                        } else {
//                            // Snap end of wire to nearest connector
//                            int ex = mx, ey = my;
//                            for (auto& c : allConnectors) {
//                                if (abs(ex - c.x) <= 2 && abs(ey - c.y) <= 2) {
//                                    ex = c.x; ey = c.y;
//                                    break;
//                                }
//                            }
//                            wires.push_back({{wireX1, wireY1}, {ex, ey}});
//                            wireStartSelected = false;
//                            currentTool = NONE;
//                            showPreview = false;
//                        }
                    }
                    else if (currentTool == VSIN) {
                        int cx, cy; snapPointWithConnectorsFirst(mx, my, cx, cy);
                        PlacedElement e{};  // ← صفر-مقداردهی ایمن
                        e.type  = VSIN;
                        e.x     = cx;
                        e.y     = cy;
                        e.name  = pendingVAC.name.empty()? "VAC" : pendingVAC.name;
                        e.value = pendingVAC.offset;   // ← offset (DC) را حتماً ست کن
                        e.amp   = pendingVAC.amp;
                        e.freq  = pendingVAC.freq;
                        e.phase = pendingVAC.phase;    // فعلاً در ذخیره‌سازی استفاده نمی‌شود
                        placedElements.push_back(e);
                        // ⬇⬇⬇ این ۳ خط را اضافه کن
                        currentTool = NONE;
                        showPreview = false;
                        // (اختیاری) پاک‌سازی بافرها
                        pendingVAC = {};

                    }

                    else if (currentTool != NONE) {

                        // Snap element end if near existing connector
                        int cx = mx, cy = my; snapPointWithConnectorsFirst(mx, my, cx, cy);
                        // For two-ended elements except GND
                        const int PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB; // = 40

                        if (currentTool != GROUND && currentTool != WIRE) {
                            for (const auto& c : allConnectors) {
                                // left end of element at cx-25, cy
                                if (abs((cx - PIN_OFFSET) - c.x) <= 2 && abs(cy - c.y) <= 2) {
                                    cx = c.x + PIN_OFFSET; cy = c.y; break;
                                }
                                // right end at cx+25, cy
                                if (abs((cx + PIN_OFFSET) - c.x) <= 2 && abs(cy - c.y) <= 2) {
                                    cx = c.x - PIN_OFFSET; cy = c.y; break;
                                }
                            }
                        } else if (currentTool == GROUND) {
                            isGroundPlaced = true;
                            // single-ended: snap center if near connector
                            for (const auto& c : allConnectors) {
                                if (abs(cx - c.x) <= 2 && abs(cy - c.y) <= 2) {
                                    cx = c.x;
                                    cy = c.y;
                                    break;
                                }
                            }
                        }
                        // ask for value and place element at snapped position

                        double val = hasPendingValue ? pendingValue : 10.0;
                        string nm = (hasPendingName ? pendingName : string(""));

                        placedElements.push_back({ currentTool, cx, cy, nm, val });

                        hasPendingValue = false;
                        pendingValue = 0.0;
                        hasPendingName = false;
                        pendingName.clear();

                        currentTool = NONE;
                        showPreview = false;
                    }
                }
            }
            // TEXTINPUT
            if (e.type == SDL_TEXTINPUT) {
                if      (vacFocus == 0) vacNameText   += e.text.text;
                else if (vacFocus == 1) vacOffsetText += e.text.text; // ← جدید
                else if (vacFocus == 2) vacAmpText    += e.text.text;
                else if (vacFocus == 3) vacFreqText   += e.text.text;
                else if (vacFocus == 4) vacPhaseText  += e.text.text;
            }
                // KEYDOWN
            else if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_BACKSPACE) {
                    string* target =
                            (vacFocus==0 ? &vacNameText :
                             (vacFocus==1 ? &vacOffsetText :
                              (vacFocus==2 ? &vacAmpText  :
                               (vacFocus==3 ? &vacFreqText : &vacPhaseText))));
                    if (!target->empty()) target->pop_back();
                } else if (e.key.keysym.sym == SDLK_TAB) {
                    vacFocus = (vacFocus + 1) % 5; // ← 5 فیلد داریم
                }
            }


            // --- PATCH C1: handle SDL_TEXTINPUT for AC fields ---
            if (showRunDialog && runChoice == RUN_AC && e.type == SDL_TEXTINPUT) {
                const char* txt = e.text.text;
                string* target = nullptr;
                if      (acFocus == 0) target = &acStartText;
                else if (acFocus == 1) target = &acStopText;
                else if (acFocus == 2) target = &acPointsPerDecText;
                else if (acFocus == 3) target = &acProbeText;

                if (target) {
                    // اجازه به اعداد/نقطه/حروف واحد: kKmMuUnN و e/E
                    for (const char* p = txt; *p; ++p) {
                        char ch = *p;
                        if (isdigit((unsigned char)ch) || ch=='.' || ch=='e' || ch=='E' ||
                            ch=='k'||ch=='K'||ch=='m'||ch=='M'||ch=='u'||ch=='U'||ch=='n'||ch=='N') {
                            target->push_back(ch);
                        }
                    }
                }
            }

            // --- PATCH C2: handle SDL_KEYDOWN for AC fields ---
            if (showRunDialog && runChoice == RUN_AC && e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_BACKSPACE) {
                    string* target = nullptr;
                    if      (acFocus == 0) target = &acStartText;
                    else if (acFocus == 1) target = &acStopText;
                    else if (acFocus == 2) target = &acPointsPerDecText;
                    else if (acFocus == 3) target = &acProbeText;
                    if (target && !target->empty()) target->pop_back();
                } else if (e.key.keysym.sym == SDLK_TAB) {
                    acFocus = (acFocus + 1) % 4;
                } else if (e.key.keysym.sym == SDLK_ESCAPE) {
                    // اختیاری: بستن دیالوگ
                    // showRunDialog = false; SDL_StopTextInput();
                }
            }

            else if (showValueDialog) {
                if (e.type == SDL_TEXTINPUT) {
                    const char* txt = e.text.text;
                    if (dialogFocus == FIELD_NAME) {
                        // فقط حروف/عدد/زیرخط برای نام (سلیقه‌ای)
                        for (const char* p = txt; *p; ++p) {
                            char ch = *p;
                            if (isalnum((unsigned char)ch) || ch == '_') {
                                inputNameText.push_back(ch);
                            }
                        }
                    } else { // FIELD_VALUE
                        // کاراکترهای رایج مقدار: اعداد، اعشار، منفی، e/E و prefixed (k,m,u,n,M,K)
                        for (const char* p = txt; *p; ++p) {
                            char ch = *p;
                            if (isdigit((unsigned char)ch) || ch=='.' || ch=='-' || ch=='e' || ch=='E' ||
                                ch=='k'||ch=='K'||ch=='m'||ch=='u'||ch=='U'||ch=='n'||ch=='M') {
                                inputValueText.push_back(ch);
                            }
                        }
                    }
                } else if (e.type == SDL_KEYDOWN) {
                    if (e.key.keysym.sym == SDLK_TAB) {
                        dialogFocus = (dialogFocus == FIELD_NAME) ? FIELD_VALUE : FIELD_NAME;
                    } else if (e.key.keysym.sym == SDLK_BACKSPACE) {
                        string &target = (dialogFocus == FIELD_NAME) ? inputNameText : inputValueText;
                        if (!target.empty()) target.pop_back();
                    } else if (e.key.keysym.sym == SDLK_RETURN) {
                        try {
                            if (inputNameText.empty()) throw runtime_error("empty name");
                            double val = parseNumber(inputValueText);
                            pendingValue = val;
                            hasPendingValue = true;
                            pendingName = inputNameText;
                            hasPendingName = true;

                            showValueDialog = false;
                            SDL_StopTextInput();
                            currentTool = selectedToolForDialog;
                            selectedToolForDialog = NONE;
                            showPreview = true;
                        } catch (...) {
                            cout << "Invalid name/value!" << endl;
                        }
                    } else if (e.key.keysym.sym == SDLK_ESCAPE) {
                        showValueDialog = false;
                        SDL_StopTextInput();
                    }
                }
            }
            else if (e.type == SDL_KEYDOWN && !showValueDialog && !showRunDialog) {
                SDL_Keycode key = e.key.keysym.sym;

                if (key == SDLK_r) currentTool = RESISTOR;
                else if (key == SDLK_c) currentTool = CAPACITOR;
                else if (key == SDLK_l) currentTool = INDUCTOR;
                else if (key == SDLK_v) currentTool = VSOURCE;
                else if (key == SDLK_i) currentTool = CSOURCE;
                else if (key == SDLK_g) currentTool = GROUND;
                else if (key == SDLK_w) currentTool = WIRE;
                else continue;
                // محدوده‌ی ابزارهای غیر از GROUND/WIRE: 30..330
                if (currentTool != GROUND && currentTool != WIRE) {
                    selectedToolForDialog = currentTool;
                    currentTool = NONE;
                    showValueDialog = true;

                    inputNameText.clear();             // ← نام خالی
                    inputValueText.clear();            // ← مقدار خالی
                    dialogFocus = FIELD_NAME;          // ← فوکوس ابتدا روی نام
                    SDL_StartTextInput();              // ← شروع دریافت ورودی متن
                }

            }

        }


        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        drawGrid(renderer, winW, winH);

        drawMenuBar(renderer, font);
        drawToolbar(renderer, font);
        drawPlacedElements(renderer);
        drawWires(renderer, wires);
        if (showFileMenu) drawFileDropdown(renderer, font);
        if (showPreview) drawPreviewElement(renderer, preX, preY);
        if (showValueDialog) {
            // پس‌زمینه دیالوگ
            SDL_SetRenderDrawColor(renderer, 200, 200, 200, 255);
            SDL_RenderFillRect(renderer, &valueDialogRect);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderDrawRect(renderer, &valueDialogRect);

            SDL_Color black = {0, 0, 0};

            // ===== Labels =====
            SDL_Surface* nameLbl = TTF_RenderText_Solid(font, "Name:", black);
            if (nameLbl) {
                SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, nameLbl);
                if (tex) {
                    SDL_Rect r = {valueDialogRect.x + 10, valueDialogRect.y + 12, nameLbl->w, nameLbl->h};
                    SDL_RenderCopy(renderer, tex, nullptr, &r);
                    SDL_DestroyTexture(tex);
                }
                SDL_FreeSurface(nameLbl);
            }

            SDL_Surface* valLbl = TTF_RenderText_Solid(font, "Value:", black);
            if (valLbl) {
                SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, valLbl);
                if (tex) {
                    SDL_Rect r = {valueDialogRect.x + 10, valueDialogRect.y + 42, valLbl->w, valLbl->h};
                    SDL_RenderCopy(renderer, tex, nullptr, &r);
                    SDL_DestroyTexture(tex);
                }
                SDL_FreeSurface(valLbl);
            }

            // ===== Name box =====
            SDL_Rect nameBox = {valueDialogRect.x + 60, valueDialogRect.y + 8, valueDialogRect.w - 70, 28};
            // پرکردن داخل
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderFillRect(renderer, &nameBox);
            // کادر (پررنگ‌تر اگر فوکوس)
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderDrawRect(renderer, &nameBox);
            if (dialogFocus == FIELD_NAME) SDL_RenderDrawRect(renderer, &nameBox);

            const char* nameToShow = inputNameText.empty() ? "e.g. R1" : inputNameText.c_str();
            SDL_Surface* nameSurf = TTF_RenderText_Solid(font, nameToShow, black);
            if (nameSurf) {
                SDL_Texture* nameTex = SDL_CreateTextureFromSurface(renderer, nameSurf);
                if (nameTex) {
                    SDL_Rect r = {nameBox.x + 6, nameBox.y + 6, nameSurf->w, nameSurf->h};
                    SDL_RenderCopy(renderer, nameTex, nullptr, &r);
                    SDL_DestroyTexture(nameTex);
                }
                SDL_FreeSurface(nameSurf);
            }

            // ===== Value box =====
            SDL_Rect valBox = {valueDialogRect.x + 60, valueDialogRect.y + 38, valueDialogRect.w - 70, 28};
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderFillRect(renderer, &valBox);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderDrawRect(renderer, &valBox);
            if (dialogFocus == FIELD_VALUE) SDL_RenderDrawRect(renderer, &valBox);

            const char* valToShow = inputValueText.empty() ? "e.g. 1k, 10u, 0.1" : inputValueText.c_str();
            SDL_Surface* valSurf = TTF_RenderText_Solid(font, valToShow, black);
            if (valSurf) {
                SDL_Texture* valTex = SDL_CreateTextureFromSurface(renderer, valSurf);
                if (valTex) {
                    SDL_Rect r = {valBox.x + 6, valBox.y + 6, valSurf->w, valSurf->h};
                    SDL_RenderCopy(renderer, valTex, nullptr, &r);
                    SDL_DestroyTexture(valTex);
                }
                SDL_FreeSurface(valSurf);
            }

            // ===== OK button =====
            SDL_Rect okButton = {valueDialogRect.x + 100, valueDialogRect.y + 80, 100, 30};
            SDL_SetRenderDrawColor(renderer, 180, 180, 180, 255);
            SDL_RenderFillRect(renderer, &okButton);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderDrawRect(renderer, &okButton);

            SDL_Surface* okSurf = TTF_RenderText_Solid(font, "OK", black);
            if (okSurf) {
                SDL_Texture* okTex = SDL_CreateTextureFromSurface(renderer, okSurf);
                if (okTex) {
                    SDL_Rect r = {okButton.x + (okButton.w - okSurf->w)/2, okButton.y + (okButton.h - okSurf->h)/2, okSurf->w, okSurf->h};
                    SDL_RenderCopy(renderer, okTex, nullptr, &r);
                    SDL_DestroyTexture(okTex);
                }
                SDL_FreeSurface(okSurf);
            }
        }
        if (showVACDialog) {
            // بک‌گراند
            SDL_SetRenderDrawColor(renderer, 200,200,200,255);
            SDL_RenderFillRect(renderer, &vacDialogRect);
            SDL_SetRenderDrawColor(renderer, 0,0,0,255);
            SDL_RenderDrawRect(renderer, &vacDialogRect);

            auto drawLbl = [&](const char* txt, int x, int y){
                SDL_Color black{0,0,0};
                SDL_Surface* s = TTF_RenderText_Solid(font, txt, black);
                if (!s) return;
                SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, s);
                SDL_Rect dst{ x, y, s->w, s->h };
                SDL_RenderCopy(renderer, t, nullptr, &dst);
                SDL_DestroyTexture(t); SDL_FreeSurface(s);
            };
            auto drawBox = [&](const string& val, int x, int y, int w, int h, bool focused){
                SDL_Rect box{ x, y, w, h };
                SDL_SetRenderDrawColor(renderer, 255,255,255,255);
                SDL_RenderFillRect(renderer, &box);
                SDL_SetRenderDrawColor(renderer, focused? 0:80, focused? 120:80, focused? 255:80, 255);
                SDL_RenderDrawRect(renderer, &box);

                SDL_Color black{0,0,0};
                const char* txt = val.empty()? "e.g. 1k / 2.2m / 10u" : val.c_str();
                SDL_Surface* s = TTF_RenderText_Solid(font, txt, black);
                if (s) {
                    SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, s);
                    SDL_Rect dst{ x+6, y+(h - s->h)/2, s->w, s->h };
                    SDL_RenderCopy(renderer, t, nullptr, &dst);
                    SDL_DestroyTexture(t); SDL_FreeSurface(s);
                }
            };

            int xL = vacDialogRect.x + 14;
            int xR = vacDialogRect.x + 220;
            int y0 = vacDialogRect.y + 16;
            int gap = 42;
            int boxW = 170, boxH = 28;

            drawLbl("Name:", xL, y0);
            drawBox(vacNameText, xR, y0-4, boxW, boxH, vacFocus==0);

            drawLbl("DC Offset (V):", xL, y0 + gap);
            drawBox(vacOffsetText, xR, y0 + gap - 4, boxW, boxH, vacFocus==1);

            drawLbl("Amplitude (V):", xL, y0 + 2*gap);
            drawBox(vacAmpText,   xR, y0 + 2*gap - 4, boxW, boxH, vacFocus==2);

            drawLbl("Frequency (Hz):", xL, y0 + 3*gap);
            drawBox(vacFreqText,  xR, y0 + 3*gap - 4, boxW, boxH, vacFocus==3);

            drawLbl("Phase (deg):", xL, y0 + 4*gap);
            drawBox(vacPhaseText, xR, y0 + 4*gap - 4, boxW, boxH, vacFocus==4);


            SDL_Rect okBtn     = { vacDialogRect.x + vacDialogRect.w - 210, vacDialogRect.y + vacDialogRect.h - 40, 90, 28 };
            SDL_Rect cancelBtn = { vacDialogRect.x + vacDialogRect.w - 105, vacDialogRect.y + vacDialogRect.h - 40, 90, 28 };

            SDL_SetRenderDrawColor(renderer, 180,180,180,255);
            SDL_RenderFillRect(renderer, &okBtn);
            SDL_RenderFillRect(renderer, &cancelBtn);
            SDL_SetRenderDrawColor(renderer, 0,0,0,255);
            SDL_RenderDrawRect(renderer, &okBtn);
            SDL_RenderDrawRect(renderer, &cancelBtn);

            drawLbl("OK",     okBtn.x + 30, okBtn.y + 5);
            drawLbl("Cancel", cancelBtn.x + 18, cancelBtn.y + 5);
        }
// --- RENDER RUN DIALOG ---
        // --- PATCH B: Render Run dialog (add AC UI) ---
        if (showRunDialog) {
            SDL_SetRenderDrawColor(renderer, 220, 220, 220, 255);
            SDL_RenderFillRect(renderer, &runDialogRect);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderDrawRect(renderer, &runDialogRect);

            SDL_Color black = {0,0,0};
            auto drawText = [&](const string& s, int x, int y){
                SDL_Surface* surf = TTF_RenderText_Solid(font, s.c_str(), black);
                if (!surf) return;
                SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, surf);
                SDL_Rect r = {x, y, surf->w, surf->h};
                SDL_RenderCopy(renderer, tex, nullptr, &r);
                SDL_DestroyTexture(tex); SDL_FreeSurface(surf);
            };

            drawText("Select analysis type", runDialogRect.x + 16, runDialogRect.y + 8);

            SDL_Rect rbDC = {runDialogRect.x + 20,  runDialogRect.y + 40,  18, 18};
            SDL_Rect rbAC = {runDialogRect.x + 20,  runDialogRect.y + 70,  18, 18};

            // DC
            SDL_RenderDrawRect(renderer, &rbDC);
            if (runChoice == RUN_DC) SDL_RenderFillRect(renderer, &rbDC);
            drawText("DC (solve nodal voltages)", rbDC.x + 30, rbDC.y - 2);

            // AC
            SDL_RenderDrawRect(renderer, &rbAC);
            if (runChoice == RUN_AC) SDL_RenderFillRect(renderer, &rbAC);
            drawText("AC Sweep (log decade)", rbAC.x + 30, rbAC.y - 2);

            // فیلدهای AC (روی یک ردیف)
            SDL_Rect acStartBox = {runDialogRect.x + 180, runDialogRect.y + 65,  80, 24};
            SDL_Rect acStopBox  = {runDialogRect.x + 280, runDialogRect.y + 65,  80, 24};
            SDL_Rect acPtsBox   = {runDialogRect.x + 380, runDialogRect.y + 65,  80,  24};
            SDL_Rect acProbeBox   = {runDialogRect.x + 480, runDialogRect.y + 65,  80,  24};

            auto drawBox = [&](SDL_Rect box, const string& val, bool focused, const char* placeholder){
                SDL_SetRenderDrawColor(renderer, 255,255,255,255);
                SDL_RenderFillRect(renderer, &box);
                SDL_SetRenderDrawColor(renderer, 0,0,0,255);
                SDL_RenderDrawRect(renderer, &box);
                if (focused) SDL_RenderDrawRect(renderer, &box);

                const string& s = val.empty()? string(placeholder) : val;
                SDL_Surface* surf = TTF_RenderText_Solid(font, s.c_str(), black);
                if (surf) {
                    SDL_Texture* tex = SDL_CreateTextureFromSurface(renderer, surf);
                    SDL_Rect r = { box.x + 6, box.y + 4, surf->w, surf->h };
                    SDL_RenderCopy(renderer, tex, nullptr, &r);
                    SDL_DestroyTexture(tex); SDL_FreeSurface(surf);
                }
            };

            if (runChoice == RUN_AC) {
                drawText("Start (Hz)", acStartBox.x, acStartBox.y - 18);
                drawBox(acStartBox, acStartText, (acFocus==0), "e.g. 10");

                drawText("Stop (Hz)",  acStopBox.x,  acStopBox.y  - 18);
                drawBox(acStopBox,  acStopText,  (acFocus==1), "e.g. 10k");

                drawText("Pts/dec",    acPtsBox.x,   acPtsBox.y   - 18);
                drawBox(acPtsBox,   acPointsPerDecText, (acFocus==2), "e.g. 10");

                drawText("Probe Node:",    acProbeBox.x,   acProbeBox.y   - 18);
                drawBox(acProbeBox,   acProbeText, (acFocus==3), "e.g. n01");
            }

            // Buttons
            SDL_Rect okBtn  = {runDialogRect.x + runDialogRect.w - 200, runDialogRect.y + runDialogRect.h - 40, 80, 28};
            SDL_Rect cancelBtn = {okBtn.x + 100, okBtn.y, 80, 28};
            SDL_SetRenderDrawColor(renderer, 180,180,180,255);
            SDL_RenderFillRect(renderer, &okBtn);
            SDL_RenderFillRect(renderer, &cancelBtn);
            SDL_SetRenderDrawColor(renderer, 0,0,0,255);
            SDL_RenderDrawRect(renderer, &okBtn);
            SDL_RenderDrawRect(renderer, &cancelBtn);
            drawText("OK", okBtn.x + 28, okBtn.y + 5);
            drawText("Cancel", cancelBtn.x + 14, cancelBtn.y + 5);
        }





        analyzeNodeConnections();

        SDL_RenderPresent(renderer);
    }


    TTF_CloseFont(font);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    TTF_Quit();
    SDL_Quit();
    return 0;
}