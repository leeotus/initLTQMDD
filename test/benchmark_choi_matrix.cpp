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
// CSV 输出格式：
//   gate,n,choi_nodes_none,choi_nodes_sift,choi_nodes_ig,
//   ig_ratio,sift_ratio,dense_elems,dense_mb,
//   sift_ms,ig_ms,
//   pt_before_none,pt_after_none,pt_before_ig,pt_after_ig,pt_inflation_none,pt_inflation_ig

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

#include "QuantumComputation.hpp"
#include "InteractionGraph.h"
#include "algorithms/QFT.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/Entanglement.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"

using namespace std;
using Clock = chrono::high_resolution_clock;

// -----------------------------------------------------------------------
// Helper: build Choi matrix DD for a circuit
//   Λ_U = U ⊗ U*   where U* = transpose(U†)
// -----------------------------------------------------------------------
static dd::Edge buildChoi(unique_ptr<dd::Package>& pkg, dd::Edge U, unsigned short n) {
    // U* = transpose(U†)
    dd::Edge Udag = pkg->conjugateTranspose(U);
    pkg->incRef(Udag);
    dd::Edge Uconj = pkg->transpose(Udag);
    pkg->incRef(Uconj);
    pkg->decRef(Udag);
    pkg->garbageCollect();

    // Choi = kronecker(U, U*)
    // dd->kronecker treats the second argument as the "lower" subsystem
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
// Choi matrix Λ = U ⊗ U* has a natural structure:
//   - qubit i (upper) and qubit i+n (lower) are symmetric by construction
//   - two-qubit gates (a,b) in U become (a,b) + (a+n,b+n) interactions
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
    // Mirror coupling: qubit i <-> qubit i+n
    for (unsigned short q = 0; q < n; ++q) {
        unsigned short qm = (unsigned short)(q + n);
        ig.weight[q][qm]++;  ig.weight[qm][q]++;
        ig.degree[q]++;      ig.degree[qm]++;
    }
    ig.detectSymmetry();
    return ig;
}

// -----------------------------------------------------------------------
// Helper: run one reordering strategy on a DD (fresh copy via re-build)
// -----------------------------------------------------------------------
struct ReorderResult {
    unsigned int size_before;
    unsigned int size_after;
    double ms;
};

static ReorderResult reorderOnce(unique_ptr<dd::Package>& pkg, dd::Edge choi,
                                  dd::DynamicReorderingStrategy strat,
                                  unsigned short totalVars) {
    // Make a copy by re-referencing (we work on the same pkg, reorder in place)
    pkg->incRef(choi);
    qc::permutationMap vm;
    for (unsigned short q = 0; q < totalVars; ++q) vm[q] = q;

    unsigned int before = pkg->size(choi);

    auto t0 = Clock::now();
    dd::Edge reordered = pkg->dynamicReorder(choi, vm, strat);
    double ms = chrono::duration<double, milli>(Clock::now() - t0).count();

    unsigned int after = pkg->size(reordered);
    pkg->decRef(choi);
    pkg->garbageCollect();
    return {before, after, ms};
}

// -----------------------------------------------------------------------
// Helper: partial trace over half the qubits (upper n qubits of 2n-qubit system)
// Returns (nodes_before_trace, nodes_after_trace)
// -----------------------------------------------------------------------
static pair<unsigned int, unsigned int> partialTraceExperiment(
        unique_ptr<dd::Package>& pkg, dd::Edge choi, unsigned short n) {
    unsigned int before = pkg->size(choi);

    // Eliminate the upper n qubits (qubits n .. 2n-1)
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
// Experiment row
// -----------------------------------------------------------------------
struct Row {
    string gate;
    unsigned short n;
    unsigned int none_nodes;
    unsigned int sift_nodes;
    unsigned int ig_nodes;
    double sift_ms;
    double ig_ms;
    // partial trace on uncompressed Choi
    unsigned int pt_before_none;
    unsigned int pt_after_none;
    // partial trace on igGroupSifting-compressed Choi
    unsigned int pt_before_ig;
    unsigned int pt_after_ig;
};

// -----------------------------------------------------------------------
// Helper: build circuit from gate name
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
// Run one (gate, n) configuration
// -----------------------------------------------------------------------
static bool runRow(const string& gate, unsigned short n, Row& out) {
    try {
        auto qcPtr = buildCircuit(gate, n);
        if (!qcPtr) return false;

        // ---- Fresh package per row to avoid cross-contamination ----
        auto pkg = make_unique<dd::Package>();
        pkg->setMode(dd::Matrix);

        dd::Edge U = qcPtr->buildFunctionality(pkg);
        pkg->incRef(U);

        // ---- Build Choi matrix ----
        dd::Edge choi = buildChoi(pkg, U, n);
        // choi is already incRef'd inside buildChoi

        unsigned short totalVars = 2 * n;

        // ---- Experiment 1: nodes without reordering ----
        unsigned int none_nodes = pkg->size(choi);

        // ---- Experiment 2: partial trace on uncompressed Choi ----
        auto [pt_before_none, pt_after_none] = partialTraceExperiment(pkg, choi, n);

        // ---- Experiment 3a: Sifting ----
        // Rebuild choi (partialTraceExperiment may have touched ref counts)
        // We need a fresh Choi — rebuild from U
        // U was decRef'd inside buildChoi, so rebuild from scratch
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

        // ---- Experiment 3b: IGGroupSifting ----
        auto pkg3 = make_unique<dd::Package>();
        pkg3->setMode(dd::Matrix);
        dd::Edge U3 = qcPtr->buildFunctionality(pkg3);
        pkg3->incRef(U3);
        dd::Edge choi3 = buildChoi(pkg3, U3, n);

        // Build and register IG for the 2n-qubit Choi system
        dd::InteractionGraph choiIG_graph = buildChoiIG(*qcPtr, n);
        pkg3->setInteractionGraph(choiIG_graph);

        qc::permutationMap vm3;
        for (unsigned short q = 0; q < totalVars; ++q) vm3[q] = q;
        auto t2 = Clock::now();
        dd::Edge choiIG = pkg3->dynamicReorder(choi3, vm3, dd::IGGroupSifting);
        double ig_ms = chrono::duration<double, milli>(Clock::now() - t2).count();
        unsigned int ig_nodes = pkg3->size(choiIG);

        // ---- Experiment 4: partial trace on igGroupSifting-compressed Choi ----
        auto [pt_before_ig, pt_after_ig] = partialTraceExperiment(pkg3, choiIG, n);

        out = {gate, n,
               none_nodes, sift_nodes, ig_nodes,
               sift_ms, ig_ms,
               pt_before_none, pt_after_none,
               pt_before_ig, pt_after_ig};
        return true;

    } catch (const exception& e) {
        cerr << "[skip] " << gate << " n=" << n << ": " << e.what() << "\n";
        return false;
    }
}

// -----------------------------------------------------------------------
// Dense size: 4^(2n) elements × 16 bytes (complex128)
// -----------------------------------------------------------------------
static double denseMB(unsigned short n) {
    // 2^(2n) × 2^(2n) complex matrix = 4^(2n) elements × 16 bytes
    // but we store 2^(2n) × 2^(2n) = 4^(2n) entries
    double elems = 1.0;
    for (int i = 0; i < 4 * n; ++i) elems *= 4.0;  // 4^(2n) = 16^n
    return elems * 16.0 / (1024.0 * 1024.0);
}

static long long denseElems(unsigned short n) {
    // 4^(2n) entries — cap at LLONG_MAX to avoid overflow display
    if (n > 10) return -1;  // too large, just mark as N/A
    long long e = 1;
    for (int i = 0; i < 4 * n; ++i) e *= 4;
    return e;
}

// -----------------------------------------------------------------------
// End-to-end QPT pipeline:
//   Build Choi → [optional compress] → partial trace → reduced state
//
// Metrics per pipeline:
//   peak_nodes  = max(choi_nodes_after_compress, pt_result_nodes)
//   total_ms    = compress_ms + partial_trace_ms
//   final_nodes = nodes of the partial trace result
//
// CSV: gate,n,strategy,choi_nodes,pt_nodes,peak_nodes,compress_ms,pt_ms,total_ms
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

        // --- Step 1: compress (optional) ---
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

        // --- Step 2: partial trace (upper n qubits) ---
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
    unsigned short nmin = 2, nmax = 8;
    string mode = "compare";  // compare | e2e
    vector<string> gates = {"H_layer", "CNOT_chain", "QFT", "Grover", "Clifford"};

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--nmin") && i+1 < argc) nmin = (unsigned short)atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nmax") && i+1 < argc) nmax = (unsigned short)atoi(argv[++i]);
        else if (!strcmp(argv[i], "--mode") && i+1 < argc) mode = argv[++i];
    }

    // ================================================================
    // Mode: e2e — end-to-end QPT pipeline comparison
    // Task: Build Choi → compress → partial trace → reduced density matrix
    // Goal: Which strategy gives the lowest peak nodes & total time?
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
    // Mode: compare — compression ratio + partial trace inflation
    // ================================================================
    cout << "gate,n"
         << ",choi_nodes_none,choi_nodes_sift,choi_nodes_ig"
         << ",ig_ratio,sift_ratio"
         << ",dense_elems,dense_mb"
         << ",sift_ms,ig_ms"
         << ",pt_before_none,pt_after_none,pt_inflation_none"
         << ",pt_before_ig,pt_after_ig,pt_inflation_ig"
         << "\n";

    for (unsigned short n = nmin; n <= nmax; ++n) {
        for (const auto& gate : gates) {
            if (gate == "Grover" && n < 2) continue;

            Row r;
            if (!runRow(gate, n, r)) continue;

            long long de = denseElems(n);
            double dmb = denseMB(n);
            double ig_ratio  = r.none_nodes > 0 ? (double)r.ig_nodes   / r.none_nodes : 1.0;
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
                 << "," << fixed << setprecision(1) << r.sift_ms
                 << "," << fixed << setprecision(1) << r.ig_ms
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
