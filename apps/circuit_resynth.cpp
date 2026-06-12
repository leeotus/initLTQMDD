#include "QuantumComputation.hpp"
#include <chrono>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
using Cx = std::complex<double>;
using Matrix = std::vector<std::vector<Cx>>;

static constexpr double EPS = 1e-10;

// ============== Extract Matrix from DD ==============

static Matrix extractMatrix(dd::Package& dd, dd::Edge e, int nqubits) {
    int dim = 1 << nqubits;
    Matrix U(dim, vector<Cx>(dim, 0));
    for (int col = 0; col < dim; col++) {
        for (int row = 0; row < dim; row++) {
            string path(nqubits, '0');
            for (int q = 0; q < nqubits; q++) {
                int rowBit = (row >> q) & 1;
                int colBit = (col >> q) & 1;
                path[q] = '0' + (rowBit * 2 + colBit);
            }
            dd::ComplexValue val = dd.getValueByPath(e, path);
            U[row][col] = Cx(val.r, val.i);
        }
    }
    return U;
}

// ============== Givens Decomposition ==============
// Count the number of non-trivial two-level unitaries needed.
// Each two-level unitary on a 1-bit-adjacent pair maps to O(1) CX gates
// for 0-control, O(2) for 1-control, O(6) for 2-control, etc.

struct DecompStats {
    int totalGivens = 0;
    int zeroCtrl = 0;   // no control needed
    int oneCtrl = 0;    // 1 control
    int twoCtrl = 0;    // 2 controls
    int moreCtrl = 0;   // 3+ controls
    int estCX = 0;      // estimated CNOT count
    int estTotal = 0;   // estimated total gates
};

static DecompStats analyzeDecomposition(const Matrix& U, int nqubits) {
    int dim = U.size();
    Matrix work = U;
    DecompStats stats;

    for (int col = 0; col < dim - 1; col++) {
        for (int row = dim - 1; row > col; row--) {
            if (abs(work[row][col]) < EPS) continue;

            int diff = row ^ col;
            vector<int> path;
            path.push_back(row);
            int cur = row;
            for (int b = 0; b < nqubits; b++) {
                if ((diff >> b) & 1) {
                    cur ^= (1 << b);
                    path.push_back(cur);
                }
            }

            for (int step = 0; step < (int)path.size() - 1; step++) {
                int r = path[step];
                int p = path[step + 1];
                if (abs(work[r][col]) < EPS) break;

                Cx x = work[p][col];
                Cx y = work[r][col];
                double mag = sqrt(norm(x) + norm(y));
                if (mag < EPS) continue;

                Cx cosv = conj(x) / mag;
                Cx sinv = conj(y) / mag;

                for (int j = 0; j < dim; j++) {
                    Cx v0 = work[p][j];
                    Cx v1 = work[r][j];
                    work[p][j] = conj(cosv) * v0 + conj(sinv) * v1;
                    work[r][j] = -sinv * v0 + cosv * v1;
                }

                // Check if it's trivial (identity-like)
                if (abs(cosv - 1.0) < EPS && abs(sinv) < EPS) continue;

                stats.totalGivens++;
                int targetBit = __builtin_ctz(p ^ r);
                int state0 = ((p >> targetBit) & 1) == 0 ? p : r;
                int numCtrls = 0;
                for (int b = 0; b < nqubits; b++) {
                    if (b == targetBit) continue;
                    numCtrls++;
                }

                if (numCtrls == 0) { stats.zeroCtrl++; stats.estCX += 0; stats.estTotal += 1; }
                else if (numCtrls == 1) { stats.oneCtrl++; stats.estCX += 2; stats.estTotal += 5; }
                else if (numCtrls == 2) { stats.twoCtrl++; stats.estCX += 8; stats.estTotal += 17; }
                else { stats.moreCtrl++; stats.estCX += 4 * numCtrls; stats.estTotal += 8 * numCtrls; }
            }
        }
    }

    // Diagonal phases
    double phase0 = arg(work[0][0]);
    for (int i = 1; i < dim; i++) {
        double relPhase = arg(work[i][i]) - phase0;
        if (abs(relPhase) > EPS && abs(abs(relPhase) - 2*M_PI) > EPS) {
            stats.estTotal += 3; // approximate overhead per phase
            stats.estCX += 1;
        }
    }

    return stats;
}

// ============== Main ==============

static void printUsage(const char* prog) {
    fprintf(stderr, "量子电路重综合分析 (基于 DD + Sifting + Givens 分解)\n\n");
    fprintf(stderr, "用法: %s <input_circuit> [sifting_strategy]\n", prog);
    fprintf(stderr, "\n从输入电路构建 DD, 用 sifting 压缩, 分析重综合后的门数\n");
    fprintf(stderr, "\n策略: none, sifting, lb, ig, iglb, group, iggroup (默认 iggroup)\n");
    fprintf(stderr, "\n限制: <= 10 qubits (矩阵提取需 2^n x 2^n 空间)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    string inputFile = argv[1];
    string stratName = "iggroup";
    if (argc >= 3) stratName = argv[2];

    dd::DynamicReorderingStrategy strat = qc::QuantumComputation::parseDynSiftStrategy(stratName);

    cout << "=== 量子电路重综合分析 ===" << endl;
    cout << "输入: " << inputFile << endl;
    cout << "Sifting策略: " << stratName << endl;

    qc::QuantumComputation qc(inputFile);
    unsigned short nq = qc.getNqubits();
    unsigned long nops = qc.getNops();

    cout << "电路: " << nq << " qubits, " << nops << " gates" << endl;

    if (nq > 10) {
        cerr << "ERROR: 电路过大 (" << nq << " qubits). 限制 <= 10." << endl;
        return 1;
    }

    auto wallStart = chrono::high_resolution_clock::now();

    // Build DD
    auto ddptr = make_unique<dd::Package>();
    ddptr->setMode(dd::Matrix);
    dd::Edge ddEdge = qc.buildFunctionality(ddptr);
    ddptr->incRef(ddEdge);

    auto sizeOrig = ddptr->size(ddEdge);
    cout << "DD 原始: " << sizeOrig << " nodes" << endl;

    // Apply sifting
    int sizeSifted = sizeOrig;
    if (strat != dd::None && sizeOrig > 1) {
        qc::permutationMap pmap;
        for (unsigned short i = 0; i < nq; ++i) pmap[i] = i;

        dd::InteractionGraph ig;
        ig.initForNqubits(nq);
        for (const auto& op : qc) { ig.addGate(op); }

        int prev = sizeOrig;
        for (int iter = 0; iter < 5; ++iter) {
            switch (strat) {
                case dd::IGGroupSifting:
                    ddEdge = std::get<0>(ddptr->igGroupSifting(ddEdge, pmap, ig));
                    break;
                case dd::GroupSifting:
                    ddEdge = std::get<0>(ddptr->groupSifting(ddEdge, pmap, ig));
                    break;
                case dd::IGLBSifting:
                    ddEdge = std::get<0>(ddptr->igLbSifting(ddEdge, pmap, ig));
                    break;
                case dd::IGSifting:
                    ddEdge = std::get<0>(ddptr->igSifting(ddEdge, pmap, ig));
                    break;
                case dd::LBSifting:
                    ddEdge = std::get<0>(ddptr->lbSifting(ddEdge, pmap));
                    break;
                default:
                    ddEdge = std::get<0>(ddptr->sifting(ddEdge, pmap));
                    break;
            }
            int cur = ddptr->size(ddEdge);
            if (cur >= prev) break;
            prev = cur;
        }
        sizeSifted = ddptr->size(ddEdge);
        cout << "DD 压缩: " << sizeSifted << " nodes ("
             << fixed << setprecision(1) << (100.0 * (1.0 - (double)sizeSifted / sizeOrig))
             << "% 减少)" << endl;
    }

    // Extract matrix and analyze decomposition
    int dim = 1 << nq;
    cout << "提取酉矩阵 (" << dim << "x" << dim << ")..." << endl;
    Matrix U = extractMatrix(*ddptr, ddEdge, nq);

    cout << "Givens 分解分析中..." << endl;
    DecompStats stats = analyzeDecomposition(U, nq);

    auto wallEnd = chrono::high_resolution_clock::now();
    double totalTime = chrono::duration<double>(wallEnd - wallStart).count();

    // Theoretical minimum CX count for n-qubit unitary: (4^n - 3n - 1)/4 (lower bound)
    long long theoreticalMinCX = ((1LL << (2*nq)) - 3*nq - 1) / 4;

    cout << "---" << endl;
    cout << "原始电路:       " << nops << " gates" << endl;
    cout << "DD 节点数:      " << sizeOrig << " → " << sizeSifted
         << " (压缩 " << fixed << setprecision(1)
         << (100.0 * (1.0 - (double)sizeSifted / sizeOrig)) << "%)" << endl;
    cout << "Givens 旋转数:  " << stats.totalGivens << endl;
    cout << "  0-control:    " << stats.zeroCtrl << " (直接单 qubit 门)" << endl;
    cout << "  1-control:    " << stats.oneCtrl << " (CU, ~2 CX each)" << endl;
    cout << "  2-control:    " << stats.twoCtrl << " (C²U, ~8 CX each)" << endl;
    if (stats.moreCtrl > 0)
        cout << "  3+-control:   " << stats.moreCtrl << endl;
    cout << "预估 CX 数:     " << stats.estCX << endl;
    cout << "预估总门数:     " << stats.estTotal << endl;
    cout << "理论最优 CX:    " << theoreticalMinCX << " (n=" << nq << " qubit 下界)" << endl;
    cout << "总时间:         " << fixed << setprecision(3) << totalTime << "s" << endl;

    return 0;
}
