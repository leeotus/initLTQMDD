#include <stdio.h>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <time.h>

using namespace std;

struct Result {
    unsigned int size;
    double time;
};

Result runStrategy(qc::QuantumComputation& qc, dd::DynamicReorderingStrategy strat, bool needXorInit) {
    qc::permutationMap map = qc.initialLayout;
    auto dd = make_unique<dd::Package>();
    if (needXorInit) dd->xorInit(map);
    auto g = qc.buildFunctionality(dd);

    clock_t start = clock();
    int prev = dd->size(g);
    for (int i = 0; i < 10; ++i) {
        g = dd->dynamicReorder(g, map, strat);
        auto sz = (int)dd->size(g);
        if (sz == prev) break;
        prev = sz;
    }
    clock_t end = clock();

    return {dd->size(g), (double)(end - start) / CLOCKS_PER_SEC};
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "用法: benchmark <circuit_file>" << endl;
        return 1;
    }

    string fileName = argv[1];

    // 提取电路名
    string name = fileName;
    size_t pos = name.find_last_of('/');
    if (pos != string::npos) name = name.substr(pos + 1);
    pos = name.find_first_of('.');
    if (pos != string::npos) name = name.substr(0, pos);

    qc::QuantumComputation qc(fileName);

    // 初始DD大小
    auto ddInit = make_unique<dd::Package>();
    clock_t start = clock();
    auto gInit = qc.buildFunctionality(ddInit);
    clock_t end = clock();
    unsigned int initialSize = ddInit->size(gInit);
    double buildTime = (double)(end - start) / CLOCKS_PER_SEC;
    unsigned short nQubits = qc.getNqubits();

    // 运行各算法
    auto sift         = runStrategy(qc, dd::Sifting, false);
    auto lbSift       = runStrategy(qc, dd::LBSifting, false);
    auto upperLin     = runStrategy(qc, dd::upperLinearSifting, true);
    auto lbUpperLin   = runStrategy(qc, dd::lbUpperLinearSifting, true);
    auto lowerLin     = runStrategy(qc, dd::lowerLinearSifting, true);
    auto lbLowerLin   = runStrategy(qc, dd::lbLowerLinearSifting, true);
    auto mixLin       = runStrategy(qc, dd::mixLinearSifting, true);
    auto lbMixLin     = runStrategy(qc, dd::lbMixLinearSifting, true);

    // 输出CSV行 (无header)
    cout << name << ","
         << nQubits << ","
         << initialSize << ","
         << buildTime << ","
         << sift.size << "," << sift.time << ","
         << lbSift.size << "," << lbSift.time << ","
         << upperLin.size << "," << upperLin.time << ","
         << lbUpperLin.size << "," << lbUpperLin.time << ","
         << lowerLin.size << "," << lowerLin.time << ","
         << lbLowerLin.size << "," << lbLowerLin.time << ","
         << mixLin.size << "," << mixLin.time << ","
         << lbMixLin.size << "," << lbMixLin.time
         << endl;

    return 0;
}
