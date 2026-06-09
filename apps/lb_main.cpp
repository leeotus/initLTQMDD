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

    start = clock();
    auto graph = qc.buildFunctionality(ddptr);
    end = clock();

    auto initialSize = ddptr->size(graph);
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    cout << fileName << " 初始DD大小: " << initialSize << ", 构造时间: " << duration << "s" << endl;
    cout << endl;

    // ========== 普通 Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::Sifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  Sifting收敛:          " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    // ========== Lower Bound Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::LBSifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  Lower Bound Sifting收敛:       " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    // ========== upper Linear Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        dd->xorInit(map);
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::upperLinearSifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  Upper Linear Sifting收敛:   " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    // ========== Lower Linear Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        dd->xorInit(map);
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::lowerLinearSifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  Lower Linear Sifting收敛:   " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    // ========== Mix Linear Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        dd->xorInit(map);
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::mixLinearSifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  Mix Linear Sifting收敛:     " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    // ========== LB Mix Linear Sifting ==========
    {
        qc::permutationMap map = qc.initialLayout;
        auto dd = make_unique<dd::Package>();
        dd->xorInit(map);
        auto g = qc.buildFunctionality(dd);
        start = clock();
        int prev = dd->size(g);
        for (int i = 0; i < 10; ++i) {
            g = dd->dynamicReorder(g, map, dd::lbMixLinearSifting);
            auto sz = (int)dd->size(g);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        duration = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "  LB Mix Linear Sifting收敛: " << dd->size(g) << ", 时间: " << duration << "s" << endl;
    }

    cout << endl;
    cout << "=== 结果对比 ===" << endl;
    cout << "初始大小: " << initialSize << endl;

    return 0;
}
