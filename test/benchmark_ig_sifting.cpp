#include <stdio.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include "QuantumComputation.hpp"
#include "InteractionGraph.h"

using namespace std;

struct Result {
    unsigned int size;
    double time_ms;
    unsigned int exchanges;
};

Result runStrategy(qc::QuantumComputation& qc, dd::DynamicReorderingStrategy strat) {
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
        g = dd->dynamicReorder(g, map, strat);
        auto sz = (int)dd->size(g);
        if (sz == prev) break;
        prev = sz;
    }
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();
    return {dd->size(g), elapsed, dd->exchange_base_cases};
}

Result runIGLBSifting(qc::QuantumComputation& qc) {
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
        auto result = dd->igLbSifting(g, map, ig);
        g = std::get<0>(result);
        auto sz = (int)dd->size(g);
        if (sz == prev) break;
        prev = sz;
    }
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double, milli>(end - start).count();
    return {dd->size(g), elapsed, dd->exchange_base_cases};
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "用法: benchmark_ig_sifting <circuit_file>" << endl;
        return 1;
    }

    string fileName = argv[1];

    string name = fileName;
    size_t p = name.find_last_of('/');
    if (p != string::npos) name = name.substr(p + 1);
    p = name.find_first_of('.');
    if (p != string::npos) name = name.substr(0, p);

    qc::QuantumComputation qc(fileName);
    unsigned short nq = qc.getNqubits();

    auto ddInit = make_unique<dd::Package>();
    auto gInit = qc.buildFunctionality(ddInit);
    unsigned int initSize = ddInit->size(gInit);

    auto sift   = runStrategy(qc, dd::Sifting);
    auto lbSift = runStrategy(qc, dd::LBSifting);
    auto igLb   = runIGLBSifting(qc);

    // CSV: circuit,qubits,init_size, sift_size,sift_time,sift_ex, lb_size,lb_time,lb_ex, iglb_size,iglb_time,iglb_ex
    cout << name << ","
         << nq << ","
         << initSize << ","
         << sift.size << "," << fixed << setprecision(2) << sift.time_ms << "," << sift.exchanges << ","
         << lbSift.size << "," << fixed << setprecision(2) << lbSift.time_ms << "," << lbSift.exchanges << ","
         << igLb.size << "," << fixed << setprecision(2) << igLb.time_ms << "," << igLb.exchanges
         << endl;

    return 0;
}
