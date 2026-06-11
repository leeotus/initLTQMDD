#include <stdio.h>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <time.h>
#include <cstdlib>
#include <cstring>

using namespace std;

static void printUsage(const char* prog) {
    fprintf(stderr, "用法: %s <circuit_file> [sifting_strategy]\n", prog);
    fprintf(stderr, "\n可选 sifting 策略 (也可通过环境变量 LTQMDD_DYN_SIFT 设置):\n");
    fprintf(stderr, "  none       - 不做动态筛选 (默认)\n");
    fprintf(stderr, "  sifting    - 标准 Sifting\n");
    fprintf(stderr, "  lb         - Lower Bound Sifting\n");
    fprintf(stderr, "  tightlb    - Tight Lower Bound Sifting\n");
    fprintf(stderr, "  iglb       - IG + LB Sifting\n");
    fprintf(stderr, "  ig         - IG Sifting (无剪枝)\n");
    fprintf(stderr, "  upperls    - Upper Linear Sifting\n");
    fprintf(stderr, "  lowerls    - Lower Linear Sifting\n");
    fprintf(stderr, "  mixls      - Mix Linear Sifting\n");
    fprintf(stderr, "  lbupperls  - LB Upper Linear Sifting\n");
    fprintf(stderr, "  lblowerls  - LB Lower Linear Sifting\n");
    fprintf(stderr, "  lbmixls    - LB Mix Linear Sifting\n");
    fprintf(stderr, "  igupperls  - IG Upper Linear Sifting\n");
    fprintf(stderr, "  iglowerls  - IG Lower Linear Sifting\n");
    fprintf(stderr, "  igmixls    - IG Mix Linear Sifting\n");
    fprintf(stderr, "\n示例:\n");
    fprintf(stderr, "  %s circuit.real sifting\n", prog);
    fprintf(stderr, "  LTQMDD_DYN_SIFT=iglb %s circuit.real\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string fileName = argv[1];

    // Determine strategy: command line arg > env var > none
    std::string stratName = "none";
    if (argc >= 3) {
        stratName = argv[2];
    } else {
        const char* envStrat = getenv("LTQMDD_DYN_SIFT");
        if (envStrat && strlen(envStrat) > 0) {
            stratName = envStrat;
        }
    }

    dd::DynamicReorderingStrategy strat = qc::QuantumComputation::parseDynSiftStrategy(stratName);

    qc::QuantumComputation qc(fileName);
    auto ddptr = make_unique<dd::Package>();

    clock_t start, end;

    if (strat == dd::None) {
        // No dynamic sifting: build normally, then apply post-hoc sifting
        start = clock();
        auto graph = qc.buildFunctionality(ddptr);
        end = clock();
        double buildTime = (double)(end - start) / CLOCKS_PER_SEC;
        auto initialSize = ddptr->size(graph);
        cout << fileName << " 初始DD大小: " << initialSize << ", 构造时间: " << buildTime << "s" << endl;

        // Post-hoc sifting to convergence
        start = clock();
        qc::permutationMap map = qc.initialLayout;
        int prev = initialSize;
        for (int i = 0; i < 10; ++i) {
            graph = ddptr->dynamicReorder(graph, map, dd::Sifting);
            auto sz = (int)ddptr->size(graph);
            if (sz == prev) break;
            prev = sz;
        }
        end = clock();
        double siftTime = (double)(end - start) / CLOCKS_PER_SEC;
        cout << fileName << " sifting收敛: " << ddptr->size(graph) << ", 时间: " << siftTime << "s" << endl;

    } else {
        // Dynamic sifting during construction
        cout << fileName << " 动态筛选策略: " << stratName << endl;

        qc::permutationMap map = qc.initialLayout;
        start = clock();
        auto graph = qc.buildFunctionalityDynamic(ddptr, map, strat);
        end = clock();

        double totalTime = (double)(end - start) / CLOCKS_PER_SEC;
        auto finalSize = ddptr->size(graph);
        cout << fileName << " 最终DD大小: " << finalSize << ", 总时间: " << totalTime << "s" << endl;
    }

    return 0;
}
