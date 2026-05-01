#include <stdio.h>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <time.h>

using namespace std;
using namespace chrono;

int main(int argc, char **argv) {
    clock_t start, end; // 统计时间

    // TODO: 根据不同的输入参数来选择不同的动态化sifting算法
    std::string fileName = argv[1];
    qc::QuantumComputation qc(fileName);
    auto ddptr = make_unique<dd::Package>();
    
    // 统计构造DD的时间
    start = clock();
    auto graph = qc.buildFunctionality(ddptr);
    end = clock();

    // 获取原始的DD大小
    auto initialSize = ddptr->size(graph);
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << fileName << " " << "初始DD大小: " << initialSize << ", 构造时间: " << duration << "s\r\n";
    
    // 应用sifting算法直到收敛
    start = clock();
    int prev = initialSize;
    for(int i = 0; i < 10; ++i) {
        graph = ddptr->dynamicReorder(graph, qc.initialLayout, dd::Sifting);
        auto sz = ddptr->size(graph);
        if(sz == prev) {
            break;
        }
        prev = sz;
    }

    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << fileName << " " << "sifting直到收敛: " << ddptr->size(graph) << ", 花费时间: " << duration << "s\r\n";

    return 0;
}