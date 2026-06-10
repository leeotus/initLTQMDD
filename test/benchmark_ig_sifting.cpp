#include <stdio.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "InteractionGraph.h"

using namespace std;

struct Result {
    unsigned int size;
    double time_ms;
    unsigned int exchanges;
};

Result runSifting(qc::QuantumComputation& qc) {
    auto dd = make_unique<dd::Package>();
    auto g = qc.buildFunctionality(dd);

    qc::permutationMap map;
    unsigned short nq = qc.getNqubits();
    for (unsigned short q = 0; q < nq; ++q)
        map[q] = q;

    dd->exchange_base_cases = 0;
    auto start = chrono::high_resolution_clock::now();
    int prev = dd->size(g);
    for (int i = 0; i < 10; ++i) {
        g = dd->dynamicReorder(g, map, dd::Sifting);
        auto sz = (int)dd->size(g);
        if (sz == prev) break;
        prev = sz;
    }
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();
    return {dd->size(g), elapsed, dd->exchange_base_cases};
}

Result runIGSifting(qc::QuantumComputation& qc) {
    dd::InteractionGraph ig;
    ig.build(qc);

    auto dd = make_unique<dd::Package>();
    auto g = qc.buildFunctionality(dd);

    qc::permutationMap map;
    unsigned short nq = qc.getNqubits();
    for (unsigned short q = 0; q < nq; ++q)
        map[q] = q;

    dd->exchange_base_cases = 0;
    auto start = chrono::high_resolution_clock::now();
    int prev = dd->size(g);
    for (int i = 0; i < 10; ++i) {
        auto result = dd->igSifting(g, map, ig);
        g = std::get<0>(result);
        auto sz = (int)dd->size(g);
        if (sz == prev) break;
        prev = sz;
    }
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();
    return {dd->size(g), elapsed, dd->exchange_base_cases};
}

void printHeader() {
    cout << left
         << setw(16) << "Circuit"
         << setw(8)  << "Qubits"
         << setw(10) << "InitSize"
         << " | "
         << setw(10) << "Sift_Sz"
         << setw(12) << "Sift_T(ms)"
         << setw(10) << "Sift_Ex"
         << " | "
         << setw(10) << "IG_Sz"
         << setw(12) << "IG_T(ms)"
         << setw(10) << "IG_Ex"
         << " | "
         << setw(10) << "Improve%"
         << endl;
    cout << string(130, '-') << endl;
}

void printRow(const string& name, unsigned short nq, unsigned int initSize,
              Result sift, Result ig) {
    double improve = (sift.size > 0)
        ? (1.0 - (double)ig.size / (double)sift.size) * 100.0
        : 0.0;
    cout << left
         << setw(16) << name
         << setw(8) << nq
         << setw(10) << initSize
         << " | "
         << setw(10) << sift.size
         << setw(12) << fixed << setprecision(2) << sift.time_ms
         << setw(10) << sift.exchanges
         << " | "
         << setw(10) << ig.size
         << setw(12) << fixed << setprecision(2) << ig.time_ms
         << setw(10) << ig.exchanges
         << " | "
         << setw(10) << fixed << setprecision(2) << improve << "%"
         << endl;
}

void printIGInfo(const string& name, qc::QuantumComputation& qc) {
    dd::InteractionGraph ig;
    ig.build(qc);
    cout << "  [" << name << "] IG edges: ";
    int edgeCount = 0;
    for (int i = 0; i < ig.n; ++i)
        for (int j = i + 1; j < ig.n; ++j)
            if (ig.weight[i][j] > 0) edgeCount++;
    cout << edgeCount << ", max_degree=" << *std::max_element(ig.degree.begin(), ig.degree.end());
    cout << ", sift_order=[ ";
    auto order = ig.getSiftOrder();
    for (int i = 0; i < std::min((int)order.size(), 8); ++i)
        cout << order[i] << " ";
    if (order.size() > 8) cout << "...";
    cout << "]" << endl;
}

int main(int argc, char** argv) {
    cout << "===== Interaction Graph Sifting vs Standard Sifting Benchmark =====" << endl;
    cout << "Comparing: Sifting (baseline) vs IGSifting (IG-driven order + gravity direction)" << endl;
    cout << endl;

    if (argc > 1) {
        printHeader();
        for (int i = 1; i < argc; ++i) {
            string fileName = argv[i];
            string name = fileName;
            size_t pos = name.find_last_of('/');
            if (pos != string::npos) name = name.substr(pos + 1);
            pos = name.find_first_of('.');
            if (pos != string::npos) name = name.substr(0, pos);

            try {
                qc::QuantumComputation qc(fileName);
                unsigned short nq = qc.getNqubits();

                auto ddInit = make_unique<dd::Package>();
                auto gInit = qc.buildFunctionality(ddInit);
                unsigned int initSize = ddInit->size(gInit);

                auto sift = runSifting(qc);
                auto ig   = runIGSifting(qc);

                printRow(name, nq, initSize, sift, ig);
                printIGInfo(name, qc);
            } catch (exception& e) {
                cerr << "Error loading " << fileName << ": " << e.what() << endl;
            }
        }
    } else {
        cout << "No circuit files provided. Using built-in QFT + Grover circuits." << endl;
        cout << endl;
        printHeader();

        // QFT circuits
        for (unsigned short nq = 4; nq <= 16; nq += 2) {
            string name = "QFT_" + to_string(nq);
            try {
                qc::QFT qft(nq);

                auto ddInit = make_unique<dd::Package>();
                auto gInit = qft.buildFunctionality(ddInit);
                unsigned int initSize = ddInit->size(gInit);

                auto sift = runSifting(qft);
                auto ig   = runIGSifting(qft);

                printRow(name, nq, initSize, sift, ig);
            } catch (exception& e) {
                cerr << "Error with " << name << ": " << e.what() << endl;
            }
        }

        cout << string(130, '-') << endl;

        // Grover circuits
        for (unsigned short nq = 3; nq <= 9; nq += 2) {
            string name = "Grover_" + to_string(nq);
            try {
                qc::Grover grover(nq);

                auto ddInit = make_unique<dd::Package>();
                auto gInit = grover.buildFunctionality(ddInit);
                unsigned int initSize = ddInit->size(gInit);

                auto sift = runSifting(grover);
                auto ig   = runIGSifting(grover);

                printRow(name, nq, initSize, sift, ig);
            } catch (exception& e) {
                cerr << "Error with " << name << ": " << e.what() << endl;
            }
        }
    }

    cout << endl;
    cout << "Legend:" << endl;
    cout << "  Sift_Sz / IG_Sz    = Final DD size after convergence" << endl;
    cout << "  Sift_T / IG_T      = Total time (ms)" << endl;
    cout << "  Sift_Ex / IG_Ex    = Number of exchange operations" << endl;
    cout << "  Improve%           = (Sift_Sz - IG_Sz) / Sift_Sz * 100" << endl;
    cout << "                       Positive = IG better, Negative = Sift better" << endl;

    return 0;
}
