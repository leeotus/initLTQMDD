#include <stdio.h>
#include "QuantumComputation.hpp"
#include <chrono>
#include <cstdlib>
#include <cstring>

using namespace std;

static void printUsage(const char* prog) {
    fprintf(stderr, "量子电路等价性验证 (基于 DD + Sifting)\n\n");
    fprintf(stderr, "用法: %s <circuit_A> <circuit_B> [strategy]\n", prog);
    fprintf(stderr, "\n验证 circuit_A 与 circuit_B 是否表示相同的酉矩阵\n");
    fprintf(stderr, "原理: 构建 DD(A) 和 DD(B), 计算 A * B†, 用 sifting 压缩后检查是否为 I\n");
    fprintf(stderr, "\n可选策略 (用于 sifting 加速验证, 默认 iggroup):\n");
    fprintf(stderr, "  none       - 不做变量重排序\n");
    fprintf(stderr, "  sifting    - 标准 Sifting\n");
    fprintf(stderr, "  lb         - Lower Bound Sifting\n");
    fprintf(stderr, "  tightlb    - Tight Lower Bound Sifting\n");
    fprintf(stderr, "  ig         - IG Sifting\n");
    fprintf(stderr, "  iglb       - IG + LB Sifting\n");
    fprintf(stderr, "  group      - Group Sifting\n");
    fprintf(stderr, "  iggroup    - IG Group Sifting (推荐)\n");
    fprintf(stderr, "\n示例:\n");
    fprintf(stderr, "  %s original.real optimized.real iggroup\n", prog);
    fprintf(stderr, "  %s circA.qasm circB.qasm sifting\n", prog);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    std::string fileA = argv[1];
    std::string fileB = argv[2];

    std::string stratName = "iggroup";
    if (argc >= 4) {
        stratName = argv[3];
    }

    dd::DynamicReorderingStrategy strat = qc::QuantumComputation::parseDynSiftStrategy(stratName);

    qc::QuantumComputation qcA(fileA);
    qc::QuantumComputation qcB(fileB);

    unsigned short nqA = qcA.getNqubits();
    unsigned short nqB = qcB.getNqubits();

    if (nqA != nqB) {
        cout << "NOT EQUIVALENT: qubit count mismatch (" << nqA << " vs " << nqB << ")" << endl;
        return 2;
    }

    cout << "电路 A: " << fileA << " (" << qcA.getNops() << " ops, " << nqA << " qubits)" << endl;
    cout << "电路 B: " << fileB << " (" << qcB.getNops() << " ops, " << nqB << " qubits)" << endl;
    cout << "Sifting策略: " << stratName << endl;
    cout << "---" << endl;

    auto ddptr = make_unique<dd::Package>();
    ddptr->setMode(dd::Matrix);
    auto wallStart = std::chrono::high_resolution_clock::now();

    // Build DD(A) - no sifting, identity variable order
    dd::Edge ddA = qcA.buildFunctionality(ddptr);
    ddptr->incRef(ddA);
    auto sizeA = ddptr->size(ddA);

    auto tAfterA = std::chrono::high_resolution_clock::now();
    cout << "DD(A) 构建完成: " << sizeA << " nodes, "
         << std::chrono::duration<double>(tAfterA - wallStart).count() << "s" << endl;

    // Build DD(B) - no sifting, same package, same variable order
    dd::Edge ddB = qcB.buildFunctionality(ddptr);
    ddptr->incRef(ddB);
    auto sizeB = ddptr->size(ddB);

    auto tAfterB = std::chrono::high_resolution_clock::now();
    cout << "DD(B) 构建完成: " << sizeB << " nodes, "
         << std::chrono::duration<double>(tAfterB - tAfterA).count() << "s" << endl;

    // Compute B† and A * B†
    dd::Edge ddBdag = ddptr->conjugateTranspose(ddB);
    ddptr->incRef(ddBdag);
    ddptr->decRef(ddB);

    dd::Edge product = ddptr->multiply(ddA, ddBdag);
    ddptr->incRef(product);
    ddptr->decRef(ddA);
    ddptr->decRef(ddBdag);
    ddptr->garbageCollect();

    auto tAfterMul = std::chrono::high_resolution_clock::now();
    auto productSize = ddptr->size(product);
    cout << "A * B† 计算完成: " << productSize << " nodes, "
         << std::chrono::duration<double>(tAfterMul - tAfterB).count() << "s" << endl;

    // Apply sifting on product to compact it (this is where the strategy matters)
    if (strat != dd::None && productSize > 1) {
        qc::permutationMap productMap;
        for (unsigned short i = 0; i < nqA; ++i)
            productMap[i] = i;

        // Build combined IG for sifting guidance
        dd::InteractionGraph combinedIG;
        combinedIG.initForNqubits(nqA);
        for (const auto& op : qcA) { combinedIG.addGate(op); }
        for (const auto& op : qcB) { combinedIG.addGate(op); }

        // Apply chosen sifting strategy to convergence
        int prev = productSize;
        for (int iter = 0; iter < 5; ++iter) {
            switch (strat) {
                case dd::IGGroupSifting:
                    product = std::get<0>(ddptr->igGroupSifting(product, productMap, combinedIG));
                    break;
                case dd::GroupSifting:
                    product = std::get<0>(ddptr->groupSifting(product, productMap, combinedIG));
                    break;
                case dd::IGLBSifting:
                    product = std::get<0>(ddptr->igLbSifting(product, productMap, combinedIG));
                    break;
                case dd::IGSifting:
                    product = std::get<0>(ddptr->igSifting(product, productMap, combinedIG));
                    break;
                case dd::LBSifting:
                    product = std::get<0>(ddptr->lbSifting(product, productMap));
                    break;
                case dd::TightLBSifting:
                    product = std::get<0>(ddptr->tightLbSifting(product, productMap));
                    break;
                default:
                    product = std::get<0>(ddptr->sifting(product, productMap));
                    break;
            }
            int cur = ddptr->size(product);
            if (cur >= prev) break;
            prev = cur;
        }

        auto tAfterSift = std::chrono::high_resolution_clock::now();
        cout << "Sifting 后: " << ddptr->size(product) << " nodes, "
             << std::chrono::duration<double>(tAfterSift - tAfterMul).count() << "s" << endl;
    }

    // Check: product == identity (up to global phase)?
    dd::Edge identity = ddptr->makeIdent(0, short(nqA - 1));
    bool equivalent = false;
    if (product.p == identity.p) {
        double re = dd::ComplexNumbers::val(product.w.r);
        double im = dd::ComplexNumbers::val(product.w.i);
        double mag = re * re + im * im;
        if (std::abs(mag - 1.0) < 1e-6)
            equivalent = true;
    }

    auto wallEnd = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(wallEnd - wallStart).count();

    cout << "---" << endl;
    if (equivalent) {
        double re = dd::ComplexNumbers::val(product.w.r);
        double im = dd::ComplexNumbers::val(product.w.i);
        if (std::abs(im) < 1e-10 && std::abs(re - 1.0) < 1e-10)
            cout << "EQUIVALENT (完全等价)" << endl;
        else
            cout << "EQUIVALENT (up to global phase: " << re << " + " << im << "i)" << endl;
    } else {
        cout << "NOT EQUIVALENT" << endl;
        cout << "  A*B† DD size = " << ddptr->size(product)
             << " (identity = " << ddptr->size(identity) << ")" << endl;
    }
    cout << "总时间: " << totalTime << "s" << endl;

    return equivalent ? 0 : 2;
}
