#include <stdio.h>
#include "QuantumComputation.hpp"
#include <chrono>
#include <cstdlib>
#include <cstring>

using namespace std;

static void printUsage(const char* prog) {
    fprintf(stderr, "量子电路等价性验证 (基于 IGGroup DD)\n\n");
    fprintf(stderr, "用法: %s <circuit_A> <circuit_B> [strategy]\n", prog);
    fprintf(stderr, "\n验证 circuit_A 与 circuit_B 是否表示相同的酉矩阵\n");
    fprintf(stderr, "\n可选策略 (用于加速 DD 构建, 默认 iggroup):\n");
    fprintf(stderr, "  none       - 不做变量重排序\n");
    fprintf(stderr, "  sifting    - 标准 Sifting\n");
    fprintf(stderr, "  iglb       - IG + LB Sifting\n");
    fprintf(stderr, "  iggroup    - IG Group Sifting (推荐)\n");
    fprintf(stderr, "\n原理: 构建 U_A · U_B†, 检查是否为单位矩阵\n");
    fprintf(stderr, "\n示例:\n");
    fprintf(stderr, "  %s original.real optimized.real\n", prog);
    fprintf(stderr, "  %s circA.real circB.real iggroup\n", prog);
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

    // Load circuits
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
    cout << "策略: " << stratName << endl;
    cout << "---" << endl;

    auto ddptr = make_unique<dd::Package>();
    ddptr->setMode(dd::Matrix);
    auto wallStart = std::chrono::high_resolution_clock::now();

    // Step 1: Build DD(A) using dynamic sifting strategy
    // Then build DD(B_inverse) by reversing B's gates and taking inverses
    // Approach: Build A·B† by constructing A first, then appending B† gates
    //
    // Simpler correct approach:
    // Build A and B separately (no sifting, same package, same var order),
    // then compute product and apply sifting on the product.

    // Build DD for circuit A (no sifting, identity order)
    dd::Edge ddA = qcA.buildFunctionality(ddptr);
    ddptr->incRef(ddA);
    auto sizeA = ddptr->size(ddA);

    auto tAfterA = std::chrono::high_resolution_clock::now();
    double timeA = std::chrono::duration<double>(tAfterA - wallStart).count();
    cout << "DD(A) 构建完成: " << sizeA << " nodes, " << timeA << "s" << endl;

    // Build DD for circuit B (no sifting, same package, same var order)
    dd::Edge ddB = qcB.buildFunctionality(ddptr);
    ddptr->incRef(ddB);
    auto sizeB = ddptr->size(ddB);

    auto tAfterB = std::chrono::high_resolution_clock::now();
    double timeB = std::chrono::duration<double>(tAfterB - tAfterA).count();
    cout << "DD(B) 构建完成: " << sizeB << " nodes, " << timeB << "s" << endl;

    // Step 2: Compute B_dagger = conjugateTranspose(B)
    dd::Edge ddBdag = ddptr->conjugateTranspose(ddB);
    ddptr->incRef(ddBdag);
    ddptr->decRef(ddB);

    // Step 3: Compute product = A * B_dagger
    dd::Edge product = ddptr->multiply(ddA, ddBdag);
    ddptr->incRef(product);
    ddptr->decRef(ddA);
    ddptr->decRef(ddBdag);
    ddptr->garbageCollect();

    auto tAfterMul = std::chrono::high_resolution_clock::now();
    double timeMul = std::chrono::duration<double>(tAfterMul - tAfterB).count();
    auto productSize = ddptr->size(product);
    cout << "A * B† 计算完成: " << productSize << " nodes, " << timeMul << "s" << endl;

    // Step 4: Apply sifting on the product to compact it
    // If A==B, product is identity, sifting will reduce to minimal form quickly.
    // The IG-based sifting accelerates this compaction.
    if (strat != dd::None && productSize > 1) {
        qc::permutationMap productMap;
        for (unsigned short i = 0; i < nqA; ++i) {
            productMap[i] = i;
        }

        // Build combined IG from both circuits for sifting guidance
        dd::InteractionGraph combinedIG;
        combinedIG.initForNqubits(nqA);
        for (const auto& op : qcA) { combinedIG.addGate(op); }
        for (const auto& op : qcB) { combinedIG.addGate(op); }

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
            default:
                product = std::get<0>(ddptr->sifting(product, productMap));
                break;
        }

        auto siftedSize = ddptr->size(product);
        auto tAfterSift = std::chrono::high_resolution_clock::now();
        double timeSift = std::chrono::duration<double>(tAfterSift - tAfterMul).count();
        cout << "Sifting 后: " << siftedSize << " nodes, " << timeSift << "s" << endl;
    }

    // Step 5: Check if product is identity (up to global phase)
    // Identity matrix DD has only 1 node (the terminal) when accessed via makeIdent.
    // After normalization, identity = makeIdent with weight ONE.
    // If product.p == identity.p, they share the same DD structure.
    dd::Edge identity = ddptr->makeIdent(0, short(nqA - 1));

    bool equivalent = false;
    if (product.p == identity.p) {
        double re = dd::ComplexNumbers::val(product.w.r);
        double im = dd::ComplexNumbers::val(product.w.i);
        double mag = re * re + im * im;
        if (std::abs(mag - 1.0) < 1e-6) {
            equivalent = true;
        }
    }

    auto wallEnd = std::chrono::high_resolution_clock::now();
    double totalTime = std::chrono::duration<double>(wallEnd - wallStart).count();

    cout << "---" << endl;
    if (equivalent) {
        double re = dd::ComplexNumbers::val(product.w.r);
        double im = dd::ComplexNumbers::val(product.w.i);
        if (std::abs(im) < 1e-10 && std::abs(re - 1.0) < 1e-10) {
            cout << "EQUIVALENT (完全等价)" << endl;
        } else {
            cout << "EQUIVALENT (up to global phase: " << re << " + " << im << "i)" << endl;
        }
    } else {
        cout << "NOT EQUIVALENT" << endl;
        cout << "  A*B† DD size = " << ddptr->size(product) << " (identity = " << ddptr->size(identity) << ")" << endl;
    }
    cout << "总时间: " << totalTime << "s" << endl;

    return equivalent ? 0 : 2;
}
