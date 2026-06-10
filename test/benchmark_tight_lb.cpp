#include <stdio.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"

using namespace std;

struct Result {
    unsigned int size;
    double time_ms;
    unsigned int exchanges;
};

Result runStrategy(qc::QuantumComputation& qc, dd::DynamicReorderingStrategy strat) {
    auto dd = make_unique<dd::Package>();
    auto g = qc.buildFunctionality(dd);
    
    // 构建 identity varMap: circuit qubit i -> DD level i
    qc::permutationMap map;
    unsigned short nq = qc.getNqubits();
    for (unsigned short q = 0; q < nq; ++q) {
        map[q] = q;
    }
    dd->xorInit(map);

    dd->exchange_base_cases = 0;
    auto start = chrono::high_resolution_clock::now();
    int prev = dd->size(g);
    for (int i = 0; i < 10; ++i) {
        g = dd->dynamicReorder(g, map, strat);
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
         << setw(12) << "Circuit"
         << setw(8) << "Qubits"
         << setw(10) << "InitSize"
         << setw(12) << "Build(ms)"
         << " | "
         << setw(10) << "LB_Size"
         << setw(12) << "LB_Time(ms)"
         << setw(8) << "LB_Ex"
         << " | "
         << setw(12) << "TightLB_Sz"
         << setw(14) << "TightLB_T(ms)"
         << setw(10) << "TightLB_Ex"
         << " | "
         << setw(10) << "Sift_Size"
         << setw(12) << "Sift_T(ms)"
         << setw(8) << "Sift_Ex"
         << endl;
    cout << string(152, '-') << endl;
}

void printRow(const string& name, unsigned short nq, unsigned int initSize,
              double buildTime, Result lb, Result tight, Result sift) {
    cout << left
         << setw(12) << name
         << setw(8) << nq
         << setw(10) << initSize
         << setw(12) << fixed << setprecision(2) << buildTime
         << " | "
         << setw(10) << lb.size
         << setw(12) << fixed << setprecision(2) << lb.time_ms
         << setw(8) << lb.exchanges
         << " | "
         << setw(12) << tight.size
         << setw(14) << fixed << setprecision(2) << tight.time_ms
         << setw(10) << tight.exchanges
         << " | "
         << setw(10) << sift.size
         << setw(12) << fixed << setprecision(2) << sift.time_ms
         << setw(8) << sift.exchanges
         << endl;
}

int main(int argc, char** argv) {
    cout << "===== Tight Lower Bound vs Original Lower Bound Benchmark =====" << endl;
    cout << "Comparing: Sifting (no pruning) vs LBSifting vs TightLBSifting" << endl;
    cout << endl;

    if (argc > 1) {
        // 从文件加载电路
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
                auto buildStart = chrono::high_resolution_clock::now();
                auto gInit = qc.buildFunctionality(ddInit);
                auto buildEnd = chrono::high_resolution_clock::now();
                unsigned int initSize = ddInit->size(gInit);
                double buildTime = chrono::duration<double, milli>(buildEnd - buildStart).count();

                auto sift  = runStrategy(qc, dd::Sifting);
                auto lb    = runStrategy(qc, dd::LBSifting);
                auto tight = runStrategy(qc, dd::TightLBSifting);

                printRow(name, nq, initSize, buildTime, lb, tight, sift);
            } catch (exception& e) {
                cerr << "Error loading " << fileName << ": " << e.what() << endl;
            }
        }
    } else {
        // 内置 QFT 电路测试
        cout << "No circuit files provided. Using built-in QFT circuits." << endl;
        cout << endl;
        printHeader();

        for (unsigned short nq = 4; nq <= 14; nq += 2) {
            string name = "QFT_" + to_string(nq);
            try {
                qc::QFT qft(nq);

                auto ddInit = make_unique<dd::Package>();
                auto buildStart = chrono::high_resolution_clock::now();
                auto gInit = qft.buildFunctionality(ddInit);
                auto buildEnd = chrono::high_resolution_clock::now();
                unsigned int initSize = ddInit->size(gInit);
                double buildTime = chrono::duration<double, milli>(buildEnd - buildStart).count();

                auto sift  = runStrategy(qft, dd::Sifting);
                auto lb    = runStrategy(qft, dd::LBSifting);
                auto tight = runStrategy(qft, dd::TightLBSifting);

                printRow(name, nq, initSize, buildTime, lb, tight, sift);
            } catch (exception& e) {
                cerr << "Error with " << name << ": " << e.what() << endl;
            }
        }
    }

    cout << endl;
    cout << "Legend:" << endl;
    cout << "  LB     = LBSifting (original lower bound pruning)" << endl;
    cout << "  TightLB = TightLBSifting (tighter lower bound pruning)" << endl;
    cout << "  Sift   = Sifting (no pruning, baseline for result quality)" << endl;
    cout << "  Ex     = Number of exchange operations performed" << endl;
    cout << endl;
    cout << "Key metric: TightLB should have fewer exchanges (faster) " << endl;
    cout << "            while achieving same or similar size as LB/Sift." << endl;

    return 0;
}
