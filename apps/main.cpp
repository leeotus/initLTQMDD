#include <stdio.h>
#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/GoogleRandomCircuitSampling.hpp"
#include <time.h>
#include <chrono>
#include <cstdlib>
#include <cstring>

using namespace std;

static void printUsage(const char* prog) {
    fprintf(stderr, "用法: %s <circuit_file> [sifting_strategy] [--parallel N]\n", prog);
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
    fprintf(stderr, "  group      - Group Sifting (symmetry-aware)\n");
    fprintf(stderr, "  iggroup    - IG Group Sifting (symmetry + IG + LB)\n");
    fprintf(stderr, "\n并行构建选项:\n");
    fprintf(stderr, "  --parallel N  使用 N 个线程并行归并构建 (默认 N=4)\n");
    fprintf(stderr, "\n示例:\n");
    fprintf(stderr, "  %s circuit.real sifting\n", prog);
    fprintf(stderr, "  %s circuit.real iggroup --parallel 4\n", prog);
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
    unsigned int numThreads = 0; // 0 means not using parallel mode
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--parallel") {
            if (i + 1 < argc) {
                numThreads = std::atoi(argv[i + 1]);
                if (numThreads < 1) numThreads = 4;
                i++;
            } else {
                numThreads = 4;
            }
        } else if (stratName == "none") {
            stratName = arg;
        }
    }
    
    if (stratName == "none") {
        const char* envStrat = getenv("LTQMDD_DYN_SIFT");
        if (envStrat && strlen(envStrat) > 0) {
            stratName = envStrat;
        }
    }

    dd::DynamicReorderingStrategy strat = qc::QuantumComputation::parseDynSiftStrategy(stratName);

    qc::QuantumComputation qc(fileName);
    auto ddptr = make_unique<dd::Package>();

    if (numThreads > 0) {
        // Parallel merge build mode
        cout << fileName << " 并行归并构建, 线程数: " << numThreads << ", 最终sift策略: " << stratName << endl;

        qc::permutationMap map = qc.initialLayout;
        auto wallStart = std::chrono::high_resolution_clock::now();
        auto graph = qc.buildFunctionalityParallelMerge(ddptr, map, strat, numThreads);
        auto wallEnd = std::chrono::high_resolution_clock::now();

        double wallTime = std::chrono::duration<double>(wallEnd - wallStart).count();
        auto finalSize = ddptr->size(graph);
        cout << fileName << " [并行] 最终DD大小: " << finalSize << ", 墙钟时间: " << wallTime << "s" << endl;

    } else if (strat == dd::None) {
        // No dynamic sifting: build normally, then apply post-hoc sifting
        auto wallStart = std::chrono::high_resolution_clock::now();
        auto graph = qc.buildFunctionality(ddptr);
        auto wallEnd = std::chrono::high_resolution_clock::now();
        double buildTime = std::chrono::duration<double>(wallEnd - wallStart).count();
        auto initialSize = ddptr->size(graph);
        cout << fileName << " 初始DD大小: " << initialSize << ", 构造时间: " << buildTime << "s" << endl;

        // Post-hoc sifting to convergence
        wallStart = std::chrono::high_resolution_clock::now();
        qc::permutationMap map = qc.initialLayout;
        int prev = initialSize;
        for (int i = 0; i < 10; ++i) {
            graph = ddptr->dynamicReorder(graph, map, dd::Sifting);
            auto sz = (int)ddptr->size(graph);
            if (sz == prev) break;
            prev = sz;
        }
        wallEnd = std::chrono::high_resolution_clock::now();
        double siftTime = std::chrono::duration<double>(wallEnd - wallStart).count();
        cout << fileName << " sifting收敛: " << ddptr->size(graph) << ", 时间: " << siftTime << "s" << endl;

    } else {
        // Dynamic sifting during construction (serial)
        cout << fileName << " 动态筛选策略: " << stratName << endl;

        qc::permutationMap map = qc.initialLayout;
        auto wallStart = std::chrono::high_resolution_clock::now();
        auto graph = qc.buildFunctionalityDynamic(ddptr, map, strat);
        auto wallEnd = std::chrono::high_resolution_clock::now();

        double totalTime = std::chrono::duration<double>(wallEnd - wallStart).count();
        auto finalSize = ddptr->size(graph);
        cout << fileName << " [串行] 最终DD大小: " << finalSize << ", 墙钟时间: " << totalTime << "s" << endl;
    }

    return 0;
}
