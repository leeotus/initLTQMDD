#include <stdio.h>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <time.h>

using namespace std;
using namespace chrono;

int main(int argc, char **argv) {
    clock_t start, end;

    if (argc < 2) {
        std::cout << "用法: lb_sifting <circuit_file>" << std::endl;
        return 1;
    }

    std::string fileName = argv[1];
    qc::QuantumComputation qc(fileName);
    auto ddptr = make_unique<dd::Package>();

    // 构造DD
    start = clock();
    auto graph = qc.buildFunctionality(ddptr);
    end = clock();

    auto initialSize = ddptr->size(graph);
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << fileName << " 初始DD大小: " << initialSize << ", 构造时间: " << duration << "s" << endl;

    // 保存初始layout用于对比
    qc::permutationMap siftingMap = qc.initialLayout;
    qc::permutationMap lbSiftingMap = qc.initialLayout;

    // ========== 普通 Sifting ==========
    auto ddSift = make_unique<dd::Package>();
    auto graphSift = qc.buildFunctionality(ddSift);

    start = clock();
    int prev = ddSift->size(graphSift);
    for (int i = 0; i < 10; ++i) {
        graphSift = ddSift->dynamicReorder(graphSift, siftingMap, dd::Sifting);
        auto sz = (int)ddSift->size(graphSift);
        if (sz == prev) break;
        prev = sz;
    }
    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "  Sifting收敛: " << ddSift->size(graphSift) << ", 时间: " << duration << "s" << endl;

    // ========== Lower Bound Sifting ==========
    auto ddLB = make_unique<dd::Package>();
    auto graphLB = qc.buildFunctionality(ddLB);

    start = clock();
    prev = ddLB->size(graphLB);
    for (int i = 0; i < 10; ++i) {
        graphLB = ddLB->dynamicReorder(graphLB, lbSiftingMap, dd::LBSifting);
        auto sz = (int)ddLB->size(graphLB);
        if (sz == prev) break;
        prev = sz;
    }
    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "  LB Sifting收敛: " << ddLB->size(graphLB) << ", 时间: " << duration << "s" << endl;

    // 对比
    cout << endl;
    cout << "=== 结果对比 ===" << endl;
    cout << "初始大小:      " << initialSize << endl;
    cout << "Sifting:       " << ddSift->size(graphSift) << endl;
    cout << "LB Sifting:    " << ddLB->size(graphLB) << endl;

    return 0;
}
