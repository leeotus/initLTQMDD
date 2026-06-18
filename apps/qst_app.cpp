/**
 * QST — Quantum State Tomography
 *
 * Two modes:
 *   Mode A (Vector, N vars):
 *     |psi_true> Vector DD -> Pauli sampling -> MLE rho_recon (N-var Matrix DD)
 *     Fidelity = <psi|rho|psi>  (innerProduct after multiply)
 *   Mode B (Matrix, 2N vars):
 *     rho_true = |psi><psi|  (2N-var Matrix DD via kronecker)
 *     Projectors in 2N-var space -> MLE rho_recon (2N-var) -> F = Tr(rho_true * rho_recon)
 *
 * Usage: ./qst <circuit_file> [nMeasBases] [nIterations]
 * CSV:   Strategy,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs,Mode
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <random>
#include <map>
#include <algorithm>
#include <array>

#include "QuantumComputation.hpp"
#include "algorithms/QST.hpp"
#include "InteractionGraph.h"

using namespace std;
using namespace chrono;

// ============================================================
// Build N-qubit Pauli projector Pi(bvec, bits) as N-var Matrix DD.
// Constructs each single-qubit projector independently via makeGateDD
// (which handles identity on non-target qubits internally), then
// multiplies them together.  Using separate makeGateDD calls (each
// embedding into the full N-qubit space) avoids numerical weight
// accumulation that occurs when multiplying a running identity matrix.
// ============================================================
static dd::Edge buildProjN(unique_ptr<dd::Package>& pkg, int N,
    const vector<qc::QST::MeasurementBasis::Basis>& bvec,
    const vector<int>& bits)
{
    // Build each single-qubit projector embedded in N-qubit space
    // via makeGateDD with LINE_TARGET on qubit q only.
    // makeGateDD automatically puts identity on all other lines.
    dd::Edge proj = pkg->makeIdent(0, N - 1);
    pkg->incRef(proj);

    for (int q = 0; q < N; ++q) {
        // Get single-qubit projector matrix
        array<dd::ComplexValue, 4> mat;
        constexpr dd::ComplexValue h={0.5,0}, mh={-0.5,0}, ih={0,0.5}, mih={0,-0.5};
        auto b = bvec[q]; auto o = bits[q];
        if (b == qc::QST::MeasurementBasis::Z) {
            if (o==0){mat[0]={1,0};mat[1]={0,0};mat[2]={0,0};mat[3]={0,0};}
            else     {mat[0]={0,0};mat[1]={0,0};mat[2]={0,0};mat[3]={1,0};}
        } else if (b == qc::QST::MeasurementBasis::X) {
            if (o==0){mat[0]=h;mat[1]=h;mat[2]=h;mat[3]=h;}
            else     {mat[0]=h;mat[1]=mh;mat[2]=mh;mat[3]=h;}
        } else { // Y
            if (o==0){mat[0]=h;mat[1]=mih;mat[2]=ih;mat[3]=h;}
            else     {mat[0]=h;mat[1]=ih;mat[2]=mih;mat[3]=h;}
        }
        array<short, dd::MAXN> line; line.fill(qc::LINE_DEFAULT);
        line[q] = qc::LINE_TARGET;
        // makeGateDD with one LINE_TARGET and all others LINE_DEFAULT produces
        // exactly  I_{q-1} ⊗ mat_{q} ⊗ I_{N-q-1}  in the N-qubit space.
        dd::Edge sp = pkg->makeGateDD(mat, N, line);
        pkg->incRef(sp);
        dd::Edge tmp = pkg->multiply(proj, sp);
        pkg->incRef(tmp);
        pkg->decRef(proj); pkg->decRef(sp);
        pkg->garbageCollect();
        proj = tmp;
    }
    return proj;
}

// ============================================================
// Scale a DD by a real scalar w by modifying the root edge weight.
// Returns a new Edge with ref-count incremented (caller must decRef).
// ============================================================
static dd::Edge scaleDD(unique_ptr<dd::Package>& pkg, dd::Edge A, double w) {
    // Guard: if A is zero/null or w is invalid, return A unchanged (incRef for ownership)
    if (A.p == nullptr || A.w == dd::ComplexNumbers::ZERO || !isfinite(w) || w == 0.0) {
        pkg->incRef(A);
        return A;
    }
    // IMPORTANT: use CN::val() to correctly handle sign-bit-encoded pointers
    double a_re = dd::ComplexNumbers::val(A.w.r);
    double a_im = dd::ComplexNumbers::val(A.w.i);
    double re = a_re * w;
    double im = a_im * w;
    if (!isfinite(re)) re = 0.0;
    if (!isfinite(im)) im = 0.0;
    dd::Edge res = A;
    res.w = pkg->cn.lookup(re, im);
    pkg->incRef(res);
    return res;
}

// ============================================================
// Apply sifting
// ============================================================
// applySifting: dynamicReorder does NOT change ref counts (pass-through semantics).
// Just call it directly and update rho.
static dd::Edge applySifting(unique_ptr<dd::Package>& pkg, dd::Edge rho,
    dd::DynamicReorderingStrategy strat,
    map<unsigned short, unsigned short>& vm,
    dd::InteractionGraph& ig)
{
    if (strat == dd::None) return rho;
    int prev = pkg->size(rho);
    for (int p = 0; p < 3; ++p) {
        dd::Edge rhoS;
        if      (strat == dd::Sifting)        rhoS = pkg->dynamicReorder(rho, vm, dd::Sifting);
        else if (strat == dd::LBSifting)      rhoS = pkg->dynamicReorder(rho, vm, dd::LBSifting);
        else if (strat == dd::IGSifting)      rhoS = get<0>(pkg->igSifting(rho, vm, ig));
        else if (strat == dd::IGLBSifting)    rhoS = get<0>(pkg->igLbSifting(rho, vm, ig));
        else if (strat == dd::GroupSifting)   { auto [r,_1,_2] = pkg->groupSifting(rho,vm,ig); rhoS=r; }
        else if (strat == dd::IGGroupSifting) { auto [r,_1,_2] = pkg->igGroupSifting(rho,vm,ig); rhoS=r; }
        else rhoS = rho;
        if (rhoS.p == nullptr || rhoS.w == dd::ComplexNumbers::ZERO) break;
        rho = rhoS;
        int ns = pkg->size(rho);
        if (ns == prev) break;
        prev = ns;
    }
    return rho;
}

// ============================================================
// RhoR MLE core:
//   Given rho (nVars-var Matrix DD) and projectors {Proj_k, freq_k},
//   run nIter iterations of:
//     R = sum_k (freq_k / p_k(rho)) * Proj_k
//     rho = R * rho * R / Tr(R * rho * R)
//   Apply sifting every 3 iterations.
//   Returns final rho (incRef'd), peakDD updated.
// ============================================================
static dd::Edge runMLE(
    unique_ptr<dd::Package>& pkg,
    dd::Edge rho,
    const vector<pair<dd::Edge,double>>& projs,
    int nVars, int nIter,
    dd::DynamicReorderingStrategy strat,
    map<unsigned short,unsigned short>& vm,
    dd::InteractionGraph& ig,
    unsigned& peakDD,
    int syncInterval  = 0,  // 每隔 N 轮触发一次重排（0=关闭）
    int syncThreshold = 0)  // rho DD 节点数超过此值时触发重排（0=关闭）
                            // 两个条件满足任一即触发，均为 0 则不做同步重排
{
    for (int iter = 0; iter < nIter; ++iter) {
        // ---- 同步重排触发判断 ----
        // 条件1: 按轮次间隔（iter>0 避免首轮未初始化）
        bool byInterval  = (syncInterval  > 0) && (iter > 0) && (iter % syncInterval == 0);
        // 条件2: rho DD 节点数超过阈值
        bool byThreshold = (syncThreshold > 0) && ((int)pkg->size(rho) > syncThreshold);
        if ((byInterval || byThreshold) && strat != dd::None) {
            rho = applySifting(pkg, rho, strat, vm, ig);
        }

        // Build R
        dd::Edge R = dd::Package::DDzero;
        pkg->incRef(R);

        for (auto& [pj, freq] : projs) {
            dd::Edge rhoP = pkg->multiply(rho, pj); pkg->incRef(rhoP);
            dd::ComplexValue tv = pkg->trace(rhoP);
            pkg->decRef(rhoP);
            double pk = tv.r;

            if (!isfinite(pk) || pk < 1e-12) pk = 1e-12;
            double w  = freq / pk;
            if (!isfinite(w)) w = 1e6;
            w = min(w, 1e6);  // cap

            dd::Edge wP = scaleDD(pkg, pj, w);
            dd::Edge nR = pkg->add(R, wP);
            pkg->incRef(nR);
            pkg->decRef(R); pkg->decRef(wP);
            pkg->garbageCollect();
            R = nR;
        }

        // Guard: stop if R is degenerate
        if (R.p == nullptr || R.w == dd::ComplexNumbers::ZERO) {
            pkg->decRef(R);
            pkg->garbageCollect();
            break;
        }

        // rho' = R * rho * R
        dd::Edge rhoR   = pkg->multiply(rho, R);   pkg->incRef(rhoR);
        dd::Edge R_rhoR = pkg->multiply(R, rhoR);  pkg->incRef(R_rhoR);
        pkg->decRef(rhoR);

        // Normalize: ρ' = ρ' / Tr(ρ')
        dd::ComplexValue trv = pkg->trace(R_rhoR);
        double tr_val = trv.r;
        if (!isfinite(tr_val) || tr_val <= 0.0 || tr_val > 1e20) {
            pkg->decRef(R_rhoR);
            pkg->decRef(R);
            pkg->garbageCollect();
            break;
        }
        double inv = 1.0 / tr_val;

        dd::Edge rhoN = scaleDD(pkg, R_rhoR, inv);
        pkg->decRef(R_rhoR);

        dd::Edge oldR = rho;
        rho = rhoN;
        pkg->decRef(oldR);
        pkg->decRef(R);
        pkg->garbageCollect();

        unsigned ca = pkg->maxActive;
        if (ca > peakDD) peakDD = ca;
    }
    return rho;
}

// ============================================================
// Strategy metadata
// ============================================================
static const string   NAMES[] = {"NoReorder","Sift","LBSift","IGSift","IGLBSift","GrpSift","IGGrpSift"};
static const dd::DynamicReorderingStrategy STRATS[] = {
    dd::None, dd::Sifting, dd::LBSifting, dd::IGSifting,
    dd::IGLBSifting, dd::GroupSifting, dd::IGGroupSifting
};
static constexpr int NSTRATS = 7;

// ============================================================
// Mode A result row
// ============================================================
struct RunResult {
    double fidelity, traceDistance;
    unsigned rhoSize, peakDD;
    double timeMs;
};

// ============================================================
// Mode A: N-var QST
//   rho_recon lives in N-var Matrix DD space.
//   Projectors Pi_k are N-var.
//   p_k = Tr(rho * Pi_k)
//   Fidelity = <psi|rho|psi> = innerProduct(psi, rho*psi)
// ============================================================
static RunResult runModeA(
    unique_ptr<qc::QuantumComputation>& qc,
    const vector<pair<vector<qc::QST::MeasurementBasis::Basis>, vector<int>>>& basisSeq,
    const vector<double>& freqs,
    int N, int nIter, dd::DynamicReorderingStrategy strat,
    int syncInterval  = 0,
    int syncThreshold = 0)
{
    RunResult res{};
    auto pkg = make_unique<dd::Package>();
    pkg->setMode(dd::Matrix);  // must be Matrix mode for gate DDs and density matrices

    // Rebuild |psi_true> (Vector DD in Matrix-mode package)
    auto fM = qc->buildFunctionality(pkg); pkg->incRef(fM);
    dd::Edge i0 = pkg->makeZeroState(N);   pkg->incRef(i0);
    dd::Edge psi = pkg->multiply(fM, i0);  pkg->incRef(psi);
    pkg->decRef(i0); pkg->decRef(fM);

    // Build projectors (N-var Matrix DDs)
    vector<pair<dd::Edge,double>> projs;
    for (size_t k = 0; k < basisSeq.size(); ++k) {
        dd::Edge proj = buildProjN(pkg, N, basisSeq[k].first, basisSeq[k].second);
        projs.push_back({proj, freqs[k]});
    }

    // Initial rho = I/2^N (maximally mixed, N-var Matrix DD)
    dd::Edge rho = pkg->makeIdent(0, N - 1);
    pkg->incRef(rho);
    // Normalize: Tr(I) = 2^N, scale by 1/2^N
    {
        double inv2N = 1.0 / (double)(1 << N);
        dd::Edge rhoN = scaleDD(pkg, rho, inv2N);
        pkg->decRef(rho);
        rho = rhoN;
    }

    map<unsigned short,unsigned short> vm;
    for (int q = 0; q < N; ++q) vm[q] = (unsigned short)q;
    dd::InteractionGraph ig; ig.initForNqubits(N);
    for (int i=0;i<N;++i) for(int j=i+1;j<N;++j){ig.weight[i][j]=ig.weight[j][i]=1;ig.degree[i]++;ig.degree[j]++;}
    ig.detectSymmetry();

    unsigned peakDD = 0;
    auto t0 = high_resolution_clock::now();

    rho = runMLE(pkg, rho, projs, N, nIter, strat, vm, ig, peakDD, syncInterval, syncThreshold);

    auto t1 = high_resolution_clock::now();

    // Fidelity = <psi|rho|psi>  — compute BEFORE sifting (variable order must match psi)
    dd::Edge rhoP = pkg->multiply(rho, psi); pkg->incRef(rhoP);
    dd::ComplexValue fv = pkg->innerProduct(psi, rhoP);
    pkg->decRef(rhoP);
    double fid = max(0.0, min(1.0, fv.r));

    // Apply sifting AFTER fidelity (only affects rhoSize for benchmarking)
    if (strat != dd::None) {
        rho = applySifting(pkg, rho, strat, vm, ig);
    }

    res.fidelity      = fid;
    res.traceDistance = sqrt(max(0.0, 1.0 - fid));
    res.rhoSize       = pkg->size(rho);
    res.peakDD        = peakDD;
    res.timeMs        = duration<double,milli>(t1-t0).count();

    for (auto& [pj,_] : projs) pkg->decRef(pj);
    pkg->decRef(rho); pkg->decRef(psi);
    pkg->garbageCollect();
    return res;
}

// ============================================================
// Mode B: 2N-var QST
//   rho_true = |psi><psi| (2N-var via kronecker)
//   Projectors Π_k^{2N}: act on first N ket-vars, identity on next N bra-vars
//   p_k = Tr(rho_true * Π_k^{2N})
//   Fidelity = Tr(rho_true * rho_recon)
// ============================================================
static RunResult runModeB(
    unique_ptr<qc::QuantumComputation>& qc,
    const vector<vector<qc::QST::MeasurementBasis::Basis>>& basisList,
    int N, int nIter, dd::DynamicReorderingStrategy strat)
{
    RunResult res{};
    int N2 = N * 2;
    int nM = (int)basisList.size();

    auto pkg = make_unique<dd::Package>();
    pkg->setMode(dd::Matrix);

    auto fM = qc->buildFunctionality(pkg); pkg->incRef(fM);
    dd::Edge i0 = pkg->makeZeroState(N);   pkg->incRef(i0);
    dd::Edge psi = pkg->multiply(fM, i0);  pkg->incRef(psi);
    pkg->decRef(i0); pkg->decRef(fM);

    // rho_true = kronecker(|psi>, <psi|)
    dd::Edge bra = pkg->conjugateTranspose(psi); pkg->incRef(bra);
    dd::Edge rhoTrue = pkg->kronecker(psi, bra);  pkg->incRef(rhoTrue);
    pkg->decRef(bra);

    // For each basis, enumerate all 2^N outcomes and compute p(b,s) = Tr(rho_true * Pi_{b,s}^{2N})
    vector<pair<dd::Edge,double>> projs2N;
    for (int k = 0; k < nM; ++k) {
        const auto& bvec = basisList[k];
        int nOutcomes = 1 << N;
        for (int s = 0; s < nOutcomes; ++s) {
            vector<int> bits(N);
            for (int q = 0; q < N; ++q) bits[q] = (s >> q) & 1;

            dd::Edge proj = pkg->makeIdent(0, N2 - 1);
            pkg->incRef(proj);
            for (int q = 0; q < N; ++q) {
                array<dd::ComplexValue,4> mat;
                constexpr dd::ComplexValue h={0.5,0}, mh={-0.5,0}, ih={0,0.5}, mih={0,-0.5};
                auto b=bvec[q]; auto o=bits[q];
                if(b==qc::QST::MeasurementBasis::Z){
                    if(o==0){mat[0]={1,0};mat[1]={0,0};mat[2]={0,0};mat[3]={0,0};}
                    else    {mat[0]={0,0};mat[1]={0,0};mat[2]={0,0};mat[3]={1,0};}
                } else if(b==qc::QST::MeasurementBasis::X){
                    if(o==0){mat[0]=h;mat[1]=h;mat[2]=h;mat[3]=h;}
                    else    {mat[0]=h;mat[1]=mh;mat[2]=mh;mat[3]=h;}
                } else {
                    if(o==0){mat[0]=h;mat[1]=mih;mat[2]=ih;mat[3]=h;}
                    else    {mat[0]=h;mat[1]=ih;mat[2]=mih;mat[3]=h;}
                }
                array<short,dd::MAXN> line; line.fill(qc::LINE_DEFAULT); line[q]=qc::LINE_TARGET;
                dd::Edge sp = pkg->makeGateDD(mat, N2, line); pkg->incRef(sp);
                dd::Edge tmp = pkg->multiply(proj, sp); pkg->incRef(tmp);
                pkg->decRef(proj); pkg->decRef(sp);
                pkg->garbageCollect();
                proj = tmp;
            }

            dd::Edge rp = pkg->multiply(rhoTrue, proj); pkg->incRef(rp);
            dd::ComplexValue tv = pkg->trace(rp);
            pkg->decRef(rp);
            double prob = tv.r;
            if(prob > 1e-10) projs2N.push_back({proj, prob});
            else pkg->decRef(proj);
        }
    }

    // Normalize frequencies
    { double s=0; for(auto& [p,f]:projs2N) s+=f; if(s>0) for(auto& [p,f]:projs2N) f/=s; }

    // Initial rho = I_{2N} / 2^(2N) (maximally mixed)
    dd::Edge rho = pkg->makeIdent(0, N2 - 1);
    pkg->incRef(rho);
    {
        double inv2N2 = 1.0 / (double)(1 << N2);
        dd::Edge rhoN = scaleDD(pkg, rho, inv2N2);
        pkg->decRef(rho);
        rho = rhoN;
    }

    map<unsigned short,unsigned short> vm;
    for(int q=0;q<N2;++q) vm[q]=(unsigned short)q;
    dd::InteractionGraph ig; ig.initForNqubits(N2);
    for(int i=0;i<N2;++i) for(int j=i+1;j<N2;++j){ig.weight[i][j]=ig.weight[j][i]=1;ig.degree[i]++;ig.degree[j]++;}
    ig.detectSymmetry();

    unsigned peakDD = 0;
    auto t0 = high_resolution_clock::now();

    rho = runMLE(pkg, rho, projs2N, N2, nIter, strat, vm, ig, peakDD);

    auto t1 = high_resolution_clock::now();

    // Fidelity = Tr(rho_true * rho_recon)
    dd::Edge prod = pkg->multiply(rhoTrue, rho); pkg->incRef(prod);
    dd::ComplexValue fv = pkg->trace(prod);
    pkg->decRef(prod);
    double fid = max(0.0, min(1.0, fv.r));

    res.fidelity      = fid;
    res.traceDistance = sqrt(max(0.0, 1.0 - fid));
    res.rhoSize       = pkg->size(rho);
    res.peakDD        = peakDD;
    res.timeMs        = duration<double,milli>(t1-t0).count();

    for(auto& [pj,_]:projs2N) pkg->decRef(pj);
    pkg->decRef(rho); pkg->decRef(rhoTrue); pkg->decRef(psi);
    pkg->garbageCollect();
    return res;
}

// ============================================================
// main
// ============================================================
int main(int argc, char** argv) try {
    if (argc < 2) {
        cerr << "Usage: ./qst <circuit> [bases] [iters] [modeB=1] [syncI=0] [syncT=0] [noise=] [nprob=0.01]\n"
             << "  noise : D=depolarizing, A=amplitude damping, P=phase flip, empty=no noise\n"
             << "  nprob : noise probability per qubit per channel application\n";
        return 1;
    }

    string cf(argv[1]);
    int nM      = (argc > 2) ? atoi(argv[2]) : 10;
    int nIter   = (argc > 3) ? atoi(argv[3]) : 50;
    bool runModeB_flag = (argc > 4) ? (atoi(argv[4]) != 0) : true;
    // syncInterval:  每隔 N 轮触发同步重排（0=关闭）
    int syncInterval  = (argc > 5) ? atoi(argv[5]) : 0;
    // syncThreshold: rho DD 节点数超过此值时触发同步重排（0=关闭）
    int syncThreshold = (argc > 6) ? atoi(argv[6]) : 0;
    // noiseType: 噪声信道类型（"" = 无噪声, "D" = 去极化, "A" = 振幅阻尼, "P" = 相位翻转）
    string noiseType   = (argc > 7) ? string(argv[7]) : "";
    // noiseProb: 噪声概率（每个 qubit 的单次信道作用强度）
    double noiseProb   = (argc > 8) ? atof(argv[8]) : 0.01;

    auto qc = make_unique<qc::QuantumComputation>(cf);
    int N = (int)qc->getNqubits();

    cout << "\n=== QST: " << qc->getName() << " | N=" << N
         << " | bases=" << nM << " | iter=" << nIter
         << " | sync=(" << syncInterval << "/" << syncThreshold << ")"
         << (noiseType.empty() ? "" : (" | noise=" + noiseType + "@" + to_string(noiseProb)))
         << " ===\n\n";

    // Generate Pauli bases.
    // If nM >= 3^N, enumerate all 3^N Pauli bases (complete tomography).
    // Otherwise, randomly sample nM bases.
    mt19937 rng0(42);
    uniform_int_distribution<int> bd0(0,2);

    long long pow3N = 1;
    for(int q=0;q<N;++q) pow3N *= 3;

    vector<vector<qc::QST::MeasurementBasis::Basis>> basisList;
    if(nM >= pow3N) {
        // Enumerate all 3^N Pauli bases
        basisList.resize((size_t)pow3N);
        for(long long idx=0;idx<pow3N;++idx){
            basisList[idx].resize(N);
            long long tmp=idx;
            for(int q=0;q<N;++q){ basisList[idx][q]=static_cast<qc::QST::MeasurementBasis::Basis>(tmp%3); tmp/=3; }
        }
        cout << "  Using all " << pow3N << " Pauli bases (complete tomography)\n";
    } else {
        // Random sample
        basisList.resize(nM);
        for(int k=0;k<nM;++k){
            basisList[k].resize(N);
            for(int q=0;q<N;++q)
                basisList[k][q] = static_cast<qc::QST::MeasurementBasis::Basis>(bd0(rng0));
        }
        cout << "  Using " << nM << " random bases (partial tomography, need >= " << pow3N << " for completeness)\n";
    }
    nM = (int)basisList.size();

    // Build |psi> in reference package to compute probabilities
    auto ddRef = make_unique<dd::Package>();
    ddRef->setMode(dd::Matrix);
    auto fRef  = qc->buildFunctionality(ddRef); ddRef->incRef(fRef);
    dd::Edge i0r  = ddRef->makeZeroState(N);     ddRef->incRef(i0r);
    dd::Edge psiR = ddRef->multiply(fRef, i0r);  ddRef->incRef(psiR);
    ddRef->decRef(i0r); ddRef->decRef(fRef);

    // ── 噪声信道应用 ──
    // 含噪时 p(b,s) = Tr(Pi_{b,s} * E(|psi><psi|)) = sum_k <psi|K_k^† Pi K_k|psi>
    // 利用线性性，对每个 Kraus 算子单独计算 <psi|K^† Pi K|psi>，求和即得含噪概率。
    // 不需要显式构建密度矩阵，完全在 N 变量向量空间计算。
    // 噪声只在此处确定性地影响频率计算，MLE 迭代本身不变。
    bool usingNoise = !noiseType.empty();

    // 预构建每个 qubit 的 Kraus 算子集合（在 N 量子比特空间中嵌入）
    // noiseKrauses[q] = {K_0^(q), K_1^(q), ...}  （每个 K 是 N 变量矩阵 DD）
    vector<vector<dd::Edge>> noiseKrauses(N);
    if (usingNoise) {
        array<short,dd::MAXN> kline; kline.fill(qc::LINE_DEFAULT);
        for (int q = 0; q < N; ++q) {
            kline.fill(qc::LINE_DEFAULT); kline[q] = qc::LINE_TARGET;
            vector<array<dd::ComplexValue,4>> kMats;
            if (noiseType.find('D') != string::npos) {
                double sp = sqrt(1.0 - noiseProb), sq = sqrt(noiseProb / 3.0);
                kMats.push_back({{{sp,0},{0,0},{0,0},{sp,0}}});
                kMats.push_back({{{0,0},{sq,0},{sq,0},{0,0}}});
                kMats.push_back({{{0,0},{0,-sq},{0,sq},{0,0}}});
                kMats.push_back({{{sq,0},{0,0},{0,0},{-sq,0}}});
            } else if (noiseType.find('A') != string::npos) {
                double g = noiseProb;
                kMats.push_back({{{1,0},{0,0},{0,0},{sqrt(1-g),0}}});
                kMats.push_back({{{0,0},{sqrt(g),0},{0,0},{0,0}}});
            } else if (noiseType.find('P') != string::npos) {
                double sp = sqrt(1.0 - noiseProb), sq = sqrt(noiseProb);
                kMats.push_back({{{sp,0},{0,0},{0,0},{sp,0}}});
                kMats.push_back({{{sq,0},{0,0},{0,0},{-sq,0}}});
            }
            for (auto& m : kMats) {
                dd::Edge K = ddRef->makeGateDD(m, N, kline); ddRef->incRef(K);
                noiseKrauses[q].push_back(K);
            }
        }
        cout << "  Noise channel: type=" << noiseType << " prob=" << noiseProb << "\n";
    }

    // For each basis, enumerate all 2^N outcomes and compute p(b,s) = Tr(Pi_{b,s} * rhoTrue)
    // validSeq/freqsA accumulate (basis,outcome) pairs with p > eps
    vector<pair<vector<qc::QST::MeasurementBasis::Basis>,vector<int>>> validSeq;
    vector<double> freqsA;
    int totalOutcomes = 0;

    for(int k=0;k<nM;++k){
        int nOutcomes = 1 << N;
        totalOutcomes += nOutcomes;
        for(int s=0;s<nOutcomes;++s){
            vector<int> bits(N);
            for(int q=0;q<N;++q) bits[q] = (s >> q) & 1;

            dd::Edge proj = buildProjN(ddRef, N, basisList[k], bits);
            ddRef->incRef(proj);
            double prob = 0.0;
            if (usingNoise) {
                // 含噪: p = sum_{k1,k2,...,kN} <psi|K_{1,k1}^† ... K_{N,kN}^† Pi K_{N,kN} ... K_{1,k1}|psi>
                // 各 qubit 独立：p = sum_{k1} sum_{k2} ... sum_{kN} <psi|K^† Pi K|psi>
                // 等价于：先对 psi 施加所有 qubit 的 Kraus，再求内积
                // 简化实现：枚举所有 Kraus 组合（对少量 qubit 可行）
                // 对于 N qubit 且每 qubit m 个 Kraus：m^N 个组合
                // 此处用递归张量积方式：逐 qubit 累加
                // 构建 K_total = K_{N-1}^(N-1) ⊗ ... ⊗ K_0^(0) 并对 psi 作用
                // 简化：构建联合态 sum_{combo} K_combo|psi> 并用 innerProduct
                // 最高效实现：p = <psi| K_N^† ... K_1^† Pi K_1 ... K_N |psi>  分层累积
                vector<dd::Edge> noisyPsis; // K_combo * |psi> for all Kraus combinations
                noisyPsis.push_back(psiR); ddRef->incRef(psiR);
                // 逐 qubit 扩展 Kraus 作用
                for (int q = 0; q < N; ++q) {
                    vector<dd::Edge> newPsis;
                    for (auto& Kq : noiseKrauses[q]) {
                        for (auto& psi_prev : noisyPsis) {
                            dd::Edge psi_new = ddRef->multiply(Kq, psi_prev); ddRef->incRef(psi_new);
                            newPsis.push_back(psi_new);
                        }
                    }
                    for (auto& p : noisyPsis) ddRef->decRef(p);
                    ddRef->garbageCollect();
                    noisyPsis = newPsis;
                }
                // p = sum_combo <psi_noisy|Pi|psi_noisy> / norm
                // 注意：此处每个 K_combo|psi> 对应一个 Kraus 分支
                // 正确公式：p = sum_combo ||K_combo|psi>||^2 * <hat_psi_combo|Pi|hat_psi_combo>
                // = sum_combo <psi|K†_combo Pi K_combo|psi>
                // 用 innerProduct(noisyPsi, Pi * noisyPsi) 即可
                for (auto& psi_k : noisyPsis) {
                    dd::Edge PiPsi = ddRef->multiply(proj, psi_k); ddRef->incRef(PiPsi);
                    dd::ComplexValue ip = ddRef->innerProduct(psi_k, PiPsi);
                    ddRef->decRef(PiPsi);
                    prob += ip.r;
                }
                for (auto& p : noisyPsis) ddRef->decRef(p);
                ddRef->garbageCollect();
            } else {
                // 纯态: p = <psi|Pi|psi>
                dd::Edge Pp = ddRef->multiply(proj, psiR); ddRef->incRef(Pp);
                dd::ComplexValue ip = ddRef->innerProduct(psiR, Pp);
                ddRef->decRef(Pp);
                prob = ip.r;
            }
            ddRef->decRef(proj);
            if(prob > 1e-10){
                freqsA.push_back(prob);
                validSeq.push_back({basisList[k], bits});
            }
        }
    }
    // Normalize frequencies so sum = 1
    { double s=0; for(auto& f:freqsA) s+=f; if(s>0) for(auto& f:freqsA) f/=s; }
    ddRef->decRef(psiR);
    // 清理噪声 Kraus 算子
    for (auto& kq : noiseKrauses)
        for (auto& K : kq) ddRef->decRef(K);
    ddRef->garbageCollect();

    cout << "  Valid measurement bases: " << validSeq.size() << " / " << totalOutcomes << "\n";
    { double mn=1e9,mx=0; for(auto& f:freqsA){mn=std::min(mn,f);mx=std::max(mx,f);}
      if(!freqsA.empty()) cout << "  freq range: [" << mn << "," << mx << "]\n\n"; }

    // CSV
    string csvName=cf;
    for(auto& c:csvName) if(c=='/'||c=='\\') c='_';
    size_t dp=csvName.find_last_of('.'); if(dp!=string::npos) csvName=csvName.substr(0,dp);
    ofstream csv("qst_"+csvName+".csv");
    csv << "Strategy,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs\n";

    // Print header
    auto printHeader=[](const string& title){
        cout << "=== " << title << " ===\n";
        cout << left << setw(12) << "Strategy" << right
             << setw(10) << "Fidelity" << setw(12) << "TrDist"
             << setw(10) << "RhoSize" << setw(10) << "PeakDD"
             << setw(12) << "Time(ms)\n"
             << string(66,'-') << "\n";
    };
    auto printRow=[&](const string& name, const RunResult& r){
        cout << left << setw(12) << name << right << fixed
             << setw(10) << setprecision(4) << r.fidelity
             << setw(12) << setprecision(4) << r.traceDistance
             << setw(10) << r.rhoSize
             << setw(10) << r.peakDD
             << setw(12) << setprecision(1) << r.timeMs << "\n";
    };

    // === Mode A ===
    printHeader("Mode A: N-var rho, F=<psi|rho|psi>");
    for(int s=0;s<NSTRATS;++s){
        auto r = runModeA(qc, validSeq, freqsA, N, nIter, STRATS[s], syncInterval, syncThreshold);
        printRow(NAMES[s], r);
        csv << NAMES[s] << "," << r.fidelity << "," << r.traceDistance << ","
            << r.rhoSize << "," << r.peakDD << "," << r.timeMs << "\n";
    }
    cout << string(66,'-') << "\n\n";

    // === Mode B (only if 2N <= MAXN and flag enabled) ===
    if(runModeB_flag && N*2 <= (int)dd::MAXN) {
        printHeader("Mode B: 2N-var rho_true=|psi><psi|, F=Tr(rho_true*rho)");
        for(int s=0;s<NSTRATS;++s){
            auto r = runModeB(qc, basisList, N, nIter, STRATS[s]);
            printRow(NAMES[s], r);
        }
        cout << string(66,'-') << "\n\n";
    } else if(!runModeB_flag) {
        cout << "Mode B skipped (disabled by argument)\n\n";
    } else {
        cout << "Mode B skipped (N*2=" << N*2 << " > MAXN=" << dd::MAXN << ")\n\n";
    }

    csv.close();
    cout << "CSV: qst_" << csvName << ".csv\n";
    return 0;
} catch (const std::bad_alloc&) {
	std::cerr << "OOM: bad_alloc\r\n";
	return 42;
} catch (const std::exception &e) {
	std::cerr << "ERR: " << e.what() << "\r\n";
	return 43;
}
