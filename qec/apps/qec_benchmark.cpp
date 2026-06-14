#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <random>

#include "QECCode.hpp"
#include "SteaneCode.hpp"
#include "QuantumComputation.hpp"
#include "DDpackage.h"
#include "InteractionGraph.h"

using namespace qec;

static void printHeader() {
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << "║     QEC Benchmark - DD-based QEC Simulation Tool     ║\n";
    std::cout << "║     Steane [[7,1,3]] + QMDD + ddsim noise engine     ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n\n";
}

static std::unique_ptr<qc::QuantumComputation> buildPhysicalCircuit(
    SteaneCode& code, int numRounds)
{
    unsigned short nq = code.nPhysical() + code.nAncilla();
    auto qc = std::make_unique<qc::QuantumComputation>(nq);
    code.generateEncodingCircuitInto(*qc);
    for (int r = 0; r < numRounds; ++r)
        code.generateSyndromeExtractionInto(*qc);
    for (int i = 0; i < code.nPhysical(); ++i)
        qc->emplace_back<qc::NonUnitaryOperation>(nq, (unsigned short)i, qc::Measure);
    return qc;
}

static std::string simulateAndMeasure(
    dd::Package& dd, dd::Edge& root, unsigned short nq, std::mt19937_64& rng)
{
    std::string result(nq, '0');
    dd::Edge cur = root;
    std::uniform_real_distribution<fp> dist(0.0, 1.0);
    for (int i = nq - 1; i >= 0; --i) {
        fp p0 = dd::ComplexNumbers::mag2(cur.p->e[0].w);
        fp p1 = dd::ComplexNumbers::mag2(cur.p->e[2].w);
        fp sum = p0 + p1;
        if (sum < 1e-12) break;
        p0 /= sum;
        if (dist(rng) < p0) { cur = cur.p->e[0]; }
        else { result[cur.p->v] = '1'; cur = cur.p->e[2]; }
    }
    return result;
}

// E1: Encoding Verification
static void experiment1_encoding_verification() {
    std::cout << "══════ Experiment 1: Encoding Correctness Verification ══════\n\n";
    SteaneCode code;
    unsigned short nq = code.nPhysical() + code.nAncilla();
    auto physical = buildPhysicalCircuit(code, 0);
    std::cout << "Physical circuit: " << physical->getNops() << " gates, " << nq << " qubits\n";

    auto dd = std::make_unique<dd::Package>();
    auto root = dd->makeZeroState(nq); dd->incRef(root);
    std::array<short, qc::MAX_QUBITS> line{}; line.fill(qc::LINE_DEFAULT);
    std::map<unsigned short, unsigned short> varMap;
    for (unsigned short i = 0; i < nq; ++i) varMap[i] = i;

    for (auto& op : *physical) {
        if (op->isStandardOperation()) {
            auto dd_op = op->getDD(dd, line, varMap);
            auto tmp = dd->multiply(dd_op, root); dd->incRef(tmp); dd->decRef(root); root = tmp;
        }
    }

    std::map<std::string, unsigned int> counts;
    std::mt19937_64 rng(42);
    for (int shot = 0; shot < 1000; ++shot) {
        auto r = simulateAndMeasure(*dd, root, nq, rng); counts[r]++;
    }
    dd->decRef(root);

    int zero = 0, one = 0;
    for (auto& [s, c] : counts) {
        // Steane Z_L = Z⊗7: logical Z = XOR parity of all 7 data qubits
        int parity = 0;
        for (int i = 0; i < code.nPhysical(); ++i) if (s[i] == '1') parity ^= 1;
        if (parity == 1) one += c; else zero += c;
    }
    std::cout << "\nDecoding results (1000 shots):\n";
    std::cout << "  Logical 0: " << zero << " (" << std::fixed << std::setprecision(1) << 100.0*zero/1000.0 << "%)\n";
    std::cout << "  Logical 1: " << one << " (" << 100.0*one/1000.0 << "%)\n";
    std::cout << "  Status: " << (one == 0 ? "✓ PASSED" : "✗ FAILED") << "\n\n";
}

// E2: DD Compression
static void experiment2_dd_compression() {
    std::cout << "══════ Experiment 2: DD Compression with QEC Circuit ══════\n\n";
    SteaneCode code;
    unsigned short nq = code.nPhysical() + code.nAncilla();
    auto physical = buildPhysicalCircuit(code, 5);
    std::cout << "Circuit: " << physical->getNops() << " gates, " << nq << " qubits (5 QEC rounds)\n\n";

    std::array<short, qc::MAX_QUBITS> line{}; line.fill(qc::LINE_DEFAULT);
    std::map<unsigned short, unsigned short> varMap;
    for (unsigned short i = 0; i < nq; ++i) varMap[i] = i;

    // no sifting
    auto dd1 = std::make_unique<dd::Package>();
    auto root1 = dd1->makeZeroState(nq); dd1->incRef(root1);
    unsigned long maxNS = 0;
    for (auto& op : *physical) {
        if (op->isStandardOperation()) {
            auto dd_op = op->getDD(dd1, line, varMap);
            auto tmp = dd1->multiply(dd_op, root1); dd1->incRef(tmp); dd1->decRef(root1); root1 = tmp;
        }
        auto sz = dd1->size(root1); if (sz > maxNS) maxNS = sz;
    }
    auto sizeNS = dd1->size(root1);

    // with sifting
    auto dd2 = std::make_unique<dd::Package>();
    auto root2 = dd2->makeZeroState(nq); dd2->incRef(root2);
    unsigned long maxS = 0, reorder = 0;
    int gi = 0;
    for (auto& op : *physical) {
        if (op->isStandardOperation()) {
            if (dd2->size(root2) > 1000 && gi % 50 == 0) {
                unsigned int smin, smax;
                std::tie(root2, smin, smax) = dd2->sifting(root2, varMap); reorder++;
            }
            auto dd_op = op->getDD(dd2, line, varMap);
            auto tmp = dd2->multiply(dd_op, root2); dd2->incRef(tmp); dd2->decRef(root2); root2 = tmp;
        }
        auto sz = dd2->size(root2); if (sz > maxS) maxS = sz; gi++;
    }
    { unsigned int smin, smax; std::tie(root2, smin, smax) = dd2->sifting(root2, varMap); reorder++; }
    auto sizeS = dd2->size(root2);

    std::cout << "DD Size (no sifting):  " << sizeNS << " nodes (peak: " << maxNS << ")\n";
    std::cout << "DD Size (with sifting): " << sizeS << " nodes (peak: " << maxS << ")\n";
    std::cout << "Compression ratio:     " << std::fixed << std::setprecision(2) << (double)sizeNS/(double)sizeS << "x\n";
    std::cout << "Reordering calls:      " << reorder << "\n\n";
    dd1->decRef(root1); dd2->decRef(root2);
}

// E3: Noise Threshold
static void experiment3_noise_threshold() {
    std::cout << "══════ Experiment 3: Logical Error Rate vs Physical Error Rate ══════\n\n";
    SteaneCode code;
    unsigned short nq = code.nPhysical() + code.nAncilla();
    auto physical = buildPhysicalCircuit(code, 5);
    std::cout << "QEC rounds: 5 | Physical qubits: " << nq << " | Shots/point: 500\n\n";
    std::cout << std::left << std::setw(12) << "p_physical" << std::setw(16) << "LER" << std::setw(16) << "+/- 2σ" << std::setw(16) << "DD size" << std::setw(14) << "Time(s)\n" << std::string(74,'-') << "\n";

    std::array<short, qc::MAX_QUBITS> line{}; line.fill(qc::LINE_DEFAULT);
    std::map<unsigned short, unsigned short> varMap;
    for (unsigned short i = 0; i < nq; ++i) varMap[i] = i;

    for (int i = 0; i <= 10; ++i) {
        double p = std::pow(10.0, -5.0 + i * 0.5);
        auto t0 = std::chrono::high_resolution_clock::now();
        auto dd = std::make_unique<dd::Package>();
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<fp> ud(0.0, 1.0);
        int errors = 0; unsigned long maxSz = 0;

        for (int shot = 0; shot < 500; ++shot) {
            auto root = dd->makeZeroState(nq); dd->incRef(root);
            for (auto& op : *physical) {
                if (op->isStandardOperation()) {
                    auto dd_op = op->getDD(dd, line, varMap);
                    auto tmp = dd->multiply(dd_op, root); dd->incRef(tmp); dd->decRef(root); root = tmp;
                    if (ud(rng) < p) {
                        unsigned short tgt = op->getTargets()[0];
                        int pt = (int)(ud(rng)*3);
                        std::array<short, qc::MAX_QUBITS> nl; nl.fill(-1); nl[tgt]=2;
                        dd::Edge nop;
                        if (pt==0) nop=dd->makeGateDD(qc::Xmat,nq,nl);
                        else if (pt==1) nop=dd->makeGateDD(qc::Zmat,nq,nl);
                        else nop=dd->makeGateDD(qc::Ymat,nq,nl);
                        tmp=dd->multiply(nop,root); dd->incRef(tmp); dd->decRef(root); root=tmp;
                    }
                }
            }
            auto r = simulateAndMeasure(*dd, root, nq, rng);
            dd->decRef(root);
            int parity=0; for(int j=0;j<code.nPhysical();++j) if(r[j]=='1') parity^=1; if(parity==1) errors++;
            auto sz = dd->size(root); if (sz > maxSz) maxSz = sz;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double t = std::chrono::duration<double>(t1-t0).count();
        double ler = (double)errors/500.0;
        double s2 = 2.0*std::sqrt(ler*(1.0-ler)/500.0);
        std::cout << std::left << std::setw(12) << std::scientific << std::setprecision(2) << p
                  << std::setw(16) << std::fixed << std::setprecision(4) << ler
                  << std::setw(16) << std::setprecision(4) << s2
                  << std::setw(16) << maxSz << std::setw(14) << std::setprecision(2) << t << "\n";
    }
    std::cout << "\n";
}

// E4: Scalability
static void experiment4_round_scalability() {
    std::cout << "══════ Experiment 4: DD Size Growth vs QEC Rounds ══════\n\n";
    SteaneCode code;
    unsigned short nq = code.nPhysical() + code.nAncilla();
    std::cout << std::left << std::setw(10) << "Rounds" << std::setw(10) << "Gates" << std::setw(16) << "DD Size" << std::setw(16) << "StateVec" << std::setw(14) << "Ratio\n" << std::string(66,'-') << "\n";

    std::array<short, qc::MAX_QUBITS> line{}; line.fill(qc::LINE_DEFAULT);
    std::map<unsigned short, unsigned short> varMap;
    for (unsigned short i = 0; i < nq; ++i) varMap[i] = i;

    for (int r = 1; r <= 20; ++r) {
        auto physical = buildPhysicalCircuit(code, r);
        auto dd = std::make_unique<dd::Package>();
        auto root = dd->makeZeroState(nq); dd->incRef(root);
        for (auto& op : *physical) {
            if (op->isStandardOperation()) {
                auto dd_op = op->getDD(dd, line, varMap);
                auto tmp = dd->multiply(dd_op, root); dd->incRef(tmp); dd->decRef(root); root = tmp;
            }
        }
        auto ddSize = dd->size(root);
        auto svSize = 1ULL << nq;
        std::cout << std::left << std::setw(10) << r << std::setw(10) << physical->getNops()
                  << std::setw(16) << ddSize << std::setw(16) << svSize
                  << std::setw(14) << std::fixed << std::setprecision(2) << (double)ddSize/(double)svSize << "x\n";
        dd->decRef(root);
    }
    std::cout << "\n";
}

// E5: IG Symmetry
static void experiment5_ig_symmetry() {
    std::cout << "══════ Experiment 5: IG Symmetry Analysis for QEC Codes ══════\n\n";
    SteaneCode code;
    auto physical = buildPhysicalCircuit(code, 1);
    dd::InteractionGraph ig;
    ig.initForNqubits(code.nPhysical()+code.nAncilla());
    for (auto& op : *physical) ig.addGate(op);
    ig.detectSymmetry();

    std::cout << "Qubits: " << ig.n << " (data: 7, ancilla: 6)\n";
    std::cout << "Symmetry groups detected: " << ig.symmetricGroups.size() << "\n\n";
    for (size_t g=0; g<ig.symmetricGroups.size(); ++g) {
        std::cout << "Group " << g << " (" << ig.symmetricGroups[g].size() << " qubits): {";
        for (size_t i=0; i<ig.symmetricGroups[g].size(); ++i) {
            if (i>0) std::cout << ", ";
            int q=ig.symmetricGroups[g][i];
            std::cout << (q<7?"data":"anc") << q;
        }
        int mn=9999,mx=0;
        for (auto q:ig.symmetricGroups[g]){if(ig.degree[q]<mn)mn=ig.degree[q];if(ig.degree[q]>mx)mx=ig.degree[q];}
        std::cout << "} degree=["<<mn<<"-"<<mx<<"]\n";
    }
    std::cout << "\nIG Weight Matrix (data qubits):\n     ";
    for(int i=0;i<7;++i)std::cout<<std::setw(4)<<i; std::cout<<"\n";
    for(int i=0;i<7;++i){std::cout<<"  "<<i<<"  ";for(int j=0;j<7;++j)std::cout<<std::setw(4)<<ig.weight[i][j];std::cout<<"  deg="<<ig.degree[i]<<"\n";}
    std::cout << "\n";
}

int main() {
    printHeader();
    experiment1_encoding_verification();
    experiment2_dd_compression();
    experiment3_noise_threshold();
    experiment4_round_scalability();
    experiment5_ig_symmetry();
    std::cout << "══════ All experiments completed ══════\n";
    return 0;
}