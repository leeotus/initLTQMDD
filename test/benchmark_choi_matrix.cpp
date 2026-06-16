// benchmark_choi_matrix.cpp
// 量子过程层析（QPT）Choi 矩阵的 DD 表示 + igGroupSifting 压缩对比实验
//
// 实验维度：
//   1. Choi DD 节点数 vs Dense 理论大小（4^(2n) = 16^n 个元素）
//   2. None / Sifting / IGGroupSifting 三组压缩率对比
//   3. 多种量子门（H层、CNOT链、QFT、Grover）Choi 矩阵对比
//   4. 压缩后做 partialTrace 的中间节点膨胀量
//
// Choi 矩阵构造：
//   Λ_U = U ⊗ conj(U)   （使用 dd->kronecker(U, U_conj)，其中 U_conj = conjugateTranspose(U)^T = conj(U)）
//   实际上 Λ = (I⊗U)|Φ+><Φ+|(I⊗U†) 等价于 matrix = U ⊗ U*
//   DD 实现：先构造 U（buildFunctionality），再 kronecker(U, conjugateTranspose(U))
//   注意：conjugateTranspose 给出 U†，再取 transpose 得 U*；但更简单地直接用 kronecker(U, transpose(U†))
//
// CSV 输出格式（compare 模式）：
//   gate,n,choi_nodes_none,choi_nodes_sift,choi_nodes_ig,
//   ig_ratio,sift_ratio,dense_elems,dense_mb,
//   build_ms,sift_ms,ig_ms,pt_none_ms,pt_ig_ms,
//   pt_before_none,pt_after_none,pt_before_ig,pt_after_ig,
//   pt_inflation_none,pt_inflation_ig

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <memory>
#include <cmath>
#include <bitset>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <set>
#include <sys/resource.h>
#include <cerrno>

#include "QuantumComputation.hpp"
#include "InteractionGraph.h"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/Entanglement.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"

using namespace std;
using Clock = chrono::high_resolution_clock;

// -----------------------------------------------------------------------
// Forward declarations
// -----------------------------------------------------------------------
static pair<unsigned int, unsigned int> partialTraceExperiment(
        unique_ptr<dd::Package>& pkg, dd::Edge choi, unsigned short n);

// -----------------------------------------------------------------------
// Helper: build Choi matrix DD for a circuit
//   Λ_U = U ⊗ U*   where U* = transpose(U†)
// -----------------------------------------------------------------------
static dd::Edge buildChoi(unique_ptr<dd::Package>& pkg, dd::Edge U, unsigned short n) {
    dd::Edge Udag = pkg->conjugateTranspose(U);
    pkg->incRef(Udag);
    dd::Edge Uconj = pkg->transpose(Udag);
    pkg->incRef(Uconj);
    pkg->decRef(Udag);
    pkg->garbageCollect();

    pkg->incRef(U);
    dd::Edge choi = pkg->kronecker(U, Uconj);
    pkg->incRef(choi);

    pkg->decRef(U);
    pkg->decRef(Uconj);
    pkg->garbageCollect();
    return choi;
}

// -----------------------------------------------------------------------
// Build InteractionGraph for the Choi matrix system (2n qubits)
// -----------------------------------------------------------------------
static dd::InteractionGraph buildChoiIG(const qc::QuantumComputation& qc, unsigned short n) {
    dd::InteractionGraph ig;
    ig.initForNqubits(2 * n);

    for (const auto& op : qc) {
        vector<unsigned short> involved;
        for (auto t : op->getTargets()) involved.push_back(t);
        for (const auto& c : op->getControls()) involved.push_back(c.qubit);

        for (size_t a = 0; a < involved.size(); ++a) {
            for (size_t b = a + 1; b < involved.size(); ++b) {
                unsigned short qa = involved[a];
                unsigned short qb = involved[b];
                if (qa < n && qb < n) {
                    ig.weight[qa][qb]++;   ig.weight[qb][qa]++;
                    ig.degree[qa]++;       ig.degree[qb]++;
                    unsigned short qam = (unsigned short)(qa + n);
                    unsigned short qbm = (unsigned short)(qb + n);
                    ig.weight[qam][qbm]++; ig.weight[qbm][qam]++;
                    ig.degree[qam]++;      ig.degree[qbm]++;
                }
            }
        }
    }
    for (unsigned short q = 0; q < n; ++q) {
        unsigned short qm = (unsigned short)(q + n);
        ig.weight[q][qm]++;  ig.weight[qm][q]++;
        ig.degree[q]++;      ig.degree[qm]++;
    }
    ig.detectSymmetry();
    return ig;
}

// -----------------------------------------------------------------------
// Build circuit from gate name
// -----------------------------------------------------------------------
static unique_ptr<qc::QuantumComputation> buildCircuit(const string& gate, unsigned short n) {
    if (gate == "H_layer") {
        auto qcPtr = make_unique<qc::QuantumComputation>(n);
        for (unsigned short q = 0; q < n; ++q)
            qcPtr->emplace_back<qc::StandardOperation>(n, q, qc::H);
        return qcPtr;
    } else if (gate == "CNOT_chain") {
        auto qcPtr = make_unique<qc::QuantumComputation>(n);
        for (unsigned short q = 0; q + 1 < n; ++q)
            qcPtr->emplace_back<qc::StandardOperation>(n, qc::Control(q), (unsigned short)(q + 1), qc::X);
        return qcPtr;
    } else if (gate == "QFT") {
        return make_unique<qc::QFT>(n);
    } else if (gate == "Grover" && n >= 2) {
        return make_unique<qc::Grover>(n, 42);
    } else if (gate == "Clifford") {
        return make_unique<qc::RandomCliffordCircuit>(n, n * 3, 42);
    }
    return nullptr;
}

// -----------------------------------------------------------------------
// Partial trace over upper n qubits of a 2n-qubit system
// -----------------------------------------------------------------------
static pair<unsigned int, unsigned int> partialTraceExperiment(
        unique_ptr<dd::Package>& pkg, dd::Edge choi, unsigned short n) {
    unsigned int before = pkg->size(choi);
    bitset<dd::MAXN> elim;
    for (unsigned short q = n; q < 2 * n; ++q) elim.set(q);

    pkg->incRef(choi);
    dd::Edge reduced = pkg->partialTrace(choi, elim);
    pkg->incRef(reduced);
    unsigned int after = pkg->size(reduced);
    pkg->decRef(reduced);
    pkg->decRef(choi);
    pkg->garbageCollect();
    return {before, after};
}

// -----------------------------------------------------------------------
// Dense size
// -----------------------------------------------------------------------
static double denseMB(unsigned short n) {
    double elems = 1.0;
    for (int i = 0; i < 4 * n; ++i) elems *= 4.0;
    return elems * 16.0 / (1024.0 * 1024.0);
}

static long long denseElems(unsigned short n) {
    if (n > 10) return -1;
    long long e = 1;
    for (int i = 0; i < 4 * n; ++i) e *= 4;
    return e;
}

// -----------------------------------------------------------------------
// Run one (gate, n) configuration  (compare mode)
//   Returns true on success, fills out[none|sift|ig] columns
// -----------------------------------------------------------------------
struct CompareResult {
    string gate;
    unsigned short n;
    unsigned int none_nodes;
    unsigned int sift_nodes;
    unsigned int ig_nodes;
    double build_ms;
    double sift_ms;
    double ig_ms;
    double pt_none_ms;
    double pt_ig_ms;
    unsigned int pt_before_none;
    unsigned int pt_after_none;
    unsigned int pt_before_ig;
    unsigned int pt_after_ig;
};

static bool runCompare(const string& gate, unsigned short n, CompareResult& out) {
    try {
        auto qcPtr = buildCircuit(gate, n);
        if (!qcPtr) return false;

        unsigned short totalVars = 2 * n;

        // ---- None ----
        auto t_build = Clock::now();
        auto pkg = make_unique<dd::Package>();
        pkg->setMode(dd::Matrix);
        dd::Edge U = qcPtr->buildFunctionality(pkg);
        pkg->incRef(U);
        dd::Edge choi = buildChoi(pkg, U, n);
        double build_ms = chrono::duration<double, milli>(Clock::now() - t_build).count();

        unsigned int none_nodes = pkg->size(choi);

        // ---- partial trace on uncompressed Choi ----
        auto t_pt_none = Clock::now();
        auto [pt_before_none, pt_after_none] = partialTraceExperiment(pkg, choi, n);
        double pt_none_ms = chrono::duration<double, milli>(Clock::now() - t_pt_none).count();

        // ---- Sifting ----
        auto pkg2 = make_unique<dd::Package>();
        pkg2->setMode(dd::Matrix);
        dd::Edge U2 = qcPtr->buildFunctionality(pkg2);
        pkg2->incRef(U2);
        dd::Edge choi2 = buildChoi(pkg2, U2, n);

        qc::permutationMap vm2;
        for (unsigned short q = 0; q < totalVars; ++q) vm2[q] = q;
        auto t1 = Clock::now();
        dd::Edge choiSift = pkg2->dynamicReorder(choi2, vm2, dd::Sifting);
        double sift_ms = chrono::duration<double, milli>(Clock::now() - t1).count();
        unsigned int sift_nodes = pkg2->size(choiSift);

        // ---- IGGroupSifting ----
        auto pkg3 = make_unique<dd::Package>();
        pkg3->setMode(dd::Matrix);
        dd::Edge U3 = qcPtr->buildFunctionality(pkg3);
        pkg3->incRef(U3);
        dd::Edge choi3 = buildChoi(pkg3, U3, n);

        dd::InteractionGraph ig_graph = buildChoiIG(*qcPtr, n);
        pkg3->setInteractionGraph(ig_graph);

        qc::permutationMap vm3;
        for (unsigned short q = 0; q < totalVars; ++q) vm3[q] = q;
        auto t2 = Clock::now();
        dd::Edge choiIG = pkg3->dynamicReorder(choi3, vm3, dd::IGGroupSifting);
        double ig_ms = chrono::duration<double, milli>(Clock::now() - t2).count();
        unsigned int ig_nodes = pkg3->size(choiIG);

        // ---- partial trace on IG-compressed Choi ----
        auto t_pt_ig = Clock::now();
        auto [pt_before_ig, pt_after_ig] = partialTraceExperiment(pkg3, choiIG, n);
        double pt_ig_ms = chrono::duration<double, milli>(Clock::now() - t_pt_ig).count();

        out.gate = gate;
        out.n = n;
        out.none_nodes = none_nodes;
        out.sift_nodes = sift_nodes;
        out.ig_nodes = ig_nodes;
        out.build_ms = build_ms;
        out.sift_ms = sift_ms;
        out.ig_ms = ig_ms;
        out.pt_none_ms = pt_none_ms;
        out.pt_ig_ms = pt_ig_ms;
        out.pt_before_none = pt_before_none;
        out.pt_after_none = pt_after_none;
        out.pt_before_ig = pt_before_ig;
        out.pt_after_ig = pt_after_ig;
        return true;

    } catch (const exception& e) {
        cerr << "[skip] " << gate << " n=" << n << ": " << e.what() << "\n";
        return false;
    }
}

// -----------------------------------------------------------------------
// End-to-end QPT pipeline (e2e mode)
// -----------------------------------------------------------------------
static void runE2EPipeline(const string& gate, unsigned short n,
                           unique_ptr<qc::QuantumComputation>& qcPtr,
                           dd::DynamicReorderingStrategy strat,
                           const string& stratName) {
    try {
        auto pkg = make_unique<dd::Package>();
        pkg->setMode(dd::Matrix);

        dd::Edge U = qcPtr->buildFunctionality(pkg);
        pkg->incRef(U);
        dd::Edge choi = buildChoi(pkg, U, n);
        unsigned short totalVars = 2 * n;

        double compress_ms = 0.0;
        unsigned int choi_nodes = pkg->size(choi);

        if (strat != dd::None) {
            if (strat == dd::IGGroupSifting) {
                dd::InteractionGraph ig = buildChoiIG(*qcPtr, n);
                pkg->setInteractionGraph(ig);
            }
            qc::permutationMap vm;
            for (unsigned short q = 0; q < totalVars; ++q) vm[q] = q;
            auto t0 = Clock::now();
            choi = pkg->dynamicReorder(choi, vm, strat);
            compress_ms = chrono::duration<double, milli>(Clock::now() - t0).count();
            choi_nodes = pkg->size(choi);
        }

        bitset<dd::MAXN> elim;
        for (unsigned short q = n; q < 2 * n; ++q) elim.set(q);

        pkg->incRef(choi);
        auto t1 = Clock::now();
        dd::Edge reduced = pkg->partialTrace(choi, elim);
        double pt_ms = chrono::duration<double, milli>(Clock::now() - t1).count();
        pkg->incRef(reduced);

        unsigned int pt_nodes = pkg->size(reduced);
        unsigned int peak_nodes = max(choi_nodes, pt_nodes);
        double total_ms = compress_ms + pt_ms;

        pkg->decRef(reduced);
        pkg->decRef(choi);
        pkg->garbageCollect();

        cout << gate << "," << (int)n << "," << stratName
             << "," << choi_nodes
             << "," << pt_nodes
             << "," << peak_nodes
             << "," << fixed << setprecision(1) << compress_ms
             << "," << fixed << setprecision(1) << pt_ms
             << "," << fixed << setprecision(1) << total_ms
             << "\n";
        cout.flush();

    } catch (const exception& e) {
        cerr << "[e2e skip] " << gate << " n=" << n << " " << stratName << ": " << e.what() << "\n";
    }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv) {
    unsigned short nmin = 2, nmax = 16;
    string mode = "compare";
    vector<string> gates = {"H_layer", "CNOT_chain", "QFT", "Grover", "Clifford"};
    unsigned long memlimitMB = 0;

    // ---- Parse arguments ----
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--nmin") && i+1 < argc) {
            nmin = (unsigned short)atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--nmax") && i+1 < argc) {
            nmax = (unsigned short)atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--mode") && i+1 < argc) {
            mode = argv[++i];
        } else if (!strcmp(argv[i], "--gate") && i+1 < argc) {
            string gateStr = argv[++i];
            gates.clear();
            stringstream ss(gateStr);
            string g;
            while (getline(ss, g, ','))
                gates.push_back(g);
        } else if (!strcmp(argv[i], "--memlimit") && i+1 < argc) {
            memlimitMB = strtoul(argv[++i], nullptr, 10);
        } else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            cout << "Usage: benchmark_choi_matrix [options]\n"
                 << "Options:\n"
                 << "  --nmin N          Minimum qubit count (default: 2)\n"
                 << "  --nmax N          Maximum qubit count (default: 16)\n"
                 << "  --mode MODE       compare | e2e (default: compare)\n"
                 << "  --gate GATES      Comma-separated gate list, e.g. QFT,Grover\n"
                 << "                    Available: H_layer,CNOT_chain,QFT,Grover,Clifford\n"
                 << "  --memlimit MB     Memory limit in MB (0 = unlimited)\n"
                 << "  -h, --help        Show this help\n";
            return 0;
        }
    }

    // ---- Apply memory limit (RLIMIT_AS) ----
    if (memlimitMB > 0) {
        struct rlimit rl;
        rl.rlim_cur = rl.rlim_max = (rlim_t)memlimitMB * 1024 * 1024;
        if (setrlimit(RLIMIT_AS, &rl) != 0) {
            cerr << "[warn] setrlimit(RLIMIT_AS, " << memlimitMB << "MB) failed: "
                 << strerror(errno) << "\n";
        } else {
            cerr << "[info] Memory limit set to " << memlimitMB << " MB\n";
        }
    }

    // ================================================================
    // Mode: e2e
    // ================================================================
    if (mode == "e2e") {
        cout << "gate,n,strategy,choi_nodes,pt_nodes,peak_nodes,compress_ms,pt_ms,total_ms\n";
        for (unsigned short n = nmin; n <= nmax; ++n) {
            for (const auto& gate : gates) {
                if (gate == "Grover" && n < 2) continue;
                auto qcPtr = buildCircuit(gate, n);
                if (!qcPtr) continue;
                runE2EPipeline(gate, n, qcPtr, dd::None,           "None");
                runE2EPipeline(gate, n, qcPtr, dd::Sifting,        "Sifting");
                runE2EPipeline(gate, n, qcPtr, dd::IGGroupSifting, "IGGroupSifting");
            }
        }
        return 0;
    }

    // ================================================================
    // Mode: compare
    // ================================================================
    cout << "gate,n"
         << ",choi_nodes_none,choi_nodes_sift,choi_nodes_ig"
         << ",ig_ratio,sift_ratio"
         << ",dense_elems,dense_mb"
         << ",build_ms,sift_ms,ig_ms,pt_none_ms,pt_ig_ms"
         << ",pt_before_none,pt_after_none,pt_inflation_none"
         << ",pt_before_ig,pt_after_ig,pt_inflation_ig"
         << "\n";

    for (unsigned short n = nmin; n <= nmax; ++n) {
        for (const auto& gate : gates) {
            if (gate == "Grover" && n < 2) continue;

            CompareResult r;
            if (!runCompare(gate, n, r)) continue;

            long long de = denseElems(n);
            double dmb = denseMB(n);
            double ig_ratio   = r.none_nodes > 0 ? (double)r.ig_nodes   / r.none_nodes : 1.0;
            double sift_ratio = r.none_nodes > 0 ? (double)r.sift_nodes / r.none_nodes : 1.0;
            int pt_infl_none = (int)r.pt_after_none - (int)r.pt_before_none;
            int pt_infl_ig   = (int)r.pt_after_ig   - (int)r.pt_before_ig;

            cout << r.gate << "," << (int)r.n
                 << "," << r.none_nodes
                 << "," << r.sift_nodes
                 << "," << r.ig_nodes
                 << "," << fixed << setprecision(4) << ig_ratio
                 << "," << fixed << setprecision(4) << sift_ratio
                 << "," << (de < 0 ? -1LL : de)
                 << "," << fixed << setprecision(1) << dmb
                 << "," << fixed << setprecision(1) << r.build_ms
                 << "," << fixed << setprecision(1) << r.sift_ms
                 << "," << fixed << setprecision(1) << r.ig_ms
                 << "," << fixed << setprecision(1) << r.pt_none_ms
                 << "," << fixed << setprecision(1) << r.pt_ig_ms
                 << "," << r.pt_before_none
                 << "," << r.pt_after_none
                 << "," << pt_infl_none
                 << "," << r.pt_before_ig
                 << "," << r.pt_after_ig
                 << "," << pt_infl_ig
                 << "\n";
            cout.flush();
        }
    }
    return 0;
}
