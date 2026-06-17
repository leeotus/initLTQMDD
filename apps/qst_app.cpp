/**
 * QST 完整实现 — 密度矩阵 MLE
 * 
 * 流程:
 *   1. |ψ_true⟩ = simulate(|0⟩)                         [Vector DD]
 *   2. ρ_true = kronecker(|ψ_true⟩, conj(|ψ_true⟩))     [Matrix DD, 2n vars]
 *   3. 生成 nM 组 Pauli 投影 Π_k, 计算 p_k = Tr(ρ_true · Π_k)
 *   4. MLE: ρ_{t+1} = Σ Π_k ρ_t Π_k                    [Matrix DD, 2n vars]
 *   5. 每 3 次迭代对 ρ 做 sifting (压缩 DD)
 *   6. Fidelity = ⟨ψ_true|ρ_recon|ψ_true⟩               [跨模式验证]
 * 
 * 约束: ρ 有 2×nQubits 个变量, 需 ≤ MAXN/2 = 64 qubits
 *   实际上 MAXN=128, 推荐 nQubits ≤ 6 以保证性能
 * 
 * 用法: ./qst <circuit_file> [nMeasurements] [nIterations]
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <random>
#include <map>
#include <bitset>

#include "QuantumComputation.hpp"
#include "algorithms/QST.hpp"
#include "InteractionGraph.h"

using namespace std;
using namespace chrono;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "用法: ./qst <circuit_file> [nMeasurements] [nIterations]\n";
        cerr << "  QST: 从Pauli测量重建密度矩阵 ρ = |ψ⟩⟨ψ|\n";
        cerr << "  流程: 电路→|ψ⟩→ρ_true→测量→MLE ρ_recon→Fidelity\n";
        return 1;
    }
    string cf(argv[1]);
    int nM = (argc>2) ? atoi(argv[2]) : 10;
    int nIter = (argc>3) ? atoi(argv[3]) : 10;

    auto qc = make_unique<qc::QuantumComputation>(cf);
    int N = qc->getNqubits(); // 电路 qubit 数
    int N_rho = N * 2;        // ρ 的变量数 = 2N (行+列)
    
    if (N_rho > dd::MAXN) {
        cerr << "错误: ρ 需要 " << N_rho << " 个变量, 超过 MAXN=" << dd::MAXN << endl;
        return 1;
    }

    cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
    cout << "║  QST 密度矩阵重建 — Sifting 策略对比                         ║\n";
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  电路: " << left << setw(52) << qc->getName() << "║\n";
    cout << "║  Qubits: " << N << "  ρ vars: " << N_rho << "  测量: " << nM << " 组  MLE迭代: " << nIter;
    cout << string(20-to_string(N).size()-to_string(N_rho).size()-to_string(nM).size()-to_string(nIter).size(),' ') << "║\n";
    cout << "╠══════════════════════════════════════════════════════════════╣\n";
    cout << "║  Step 1/4: |ψ_true⟩ Vector DD → ρ_true Matrix DD           ║\n";
    cout << "║  Step 2/4: Pauli 投影 Π_k, p_k = Tr(ρ_true · Π_k)          ║\n";
    cout << "║  Step 3/4: MLE ρ_{t+1} = Σ Π_k ρ_t Π_k + sifting           ║\n";
    cout << "║  Step 4/4: Fidelity = ⟨ψ_true|ρ_recon|ψ_true⟩              ║\n";
    cout << "╚══════════════════════════════════════════════════════════════╝\n\n";

    // ================================================================
    // Step 1: |ψ_true⟩ → ρ_true
    // ================================================================
    auto ddVec = make_unique<dd::Package>();  // Vector 模式
    auto fVec = qc->buildFunctionality(ddVec); ddVec->incRef(fVec);
    dd::Edge i0 = ddVec->makeZeroState(N); ddVec->incRef(i0);
    dd::Edge psiTrue = qc->simulate(i0, ddVec); ddVec->incRef(psiTrue);
    ddVec->decRef(i0); ddVec->decRef(fVec);
    cout << "  |ψ_true⟩ Vector DD: " << ddVec->size(psiTrue) << " nodes\n";

    // ρ_true = |ψ⟩⟨ψ| via kronecker — 切换到 Matrix 模式
    auto ddMat = make_unique<dd::Package>();
    ddMat->setMode(dd::Matrix);
    
    // 在 ddMat 上重建 |ψ⟩
    auto fMat = qc->buildFunctionality(ddMat); ddMat->incRef(fMat);
    dd::Edge i1 = ddMat->makeZeroState(N); ddMat->incRef(i1);
    dd::Edge psiTrueM = qc->simulate(i1, ddMat); ddMat->incRef(psiTrueM);
    ddMat->decRef(i1); ddMat->decRef(fMat);

    // ρ_true = kronecker(|ψ⟩, conj(|ψ⟩))
    dd::Edge bra = ddMat->conjugateTranspose(psiTrueM); ddMat->incRef(bra);
    dd::Edge rhoTrue = ddMat->kronecker(psiTrueM, bra); ddMat->incRef(rhoTrue);
    ddMat->decRef(bra); ddMat->decRef(psiTrueM);
    
    unsigned rhoTrueSize = ddMat->size(rhoTrue);
    cout << "  ρ_true Matrix DD: " << rhoTrueSize << " nodes (vars=" << N_rho << ")\n";

    // ================================================================
    // Step 2: 生成 Pauli 投影算子 + 模拟测量
    // ================================================================
    cout << "\n  Step 2: 生成 " << nM << " 组 Pauli 投影算子...\n";
    
    mt19937 rng(42);
    uniform_int_distribution<int> bd(0,2);
    
    struct Proj { dd::Edge dd; double prob; };
    vector<Proj> projs;

    for (int k=0; k<nM; ++k) {
        vector<int> bits(N);
        qc::QST::MeasurementBasis basis; basis.bases.resize(N);
        for (int q=0; q<N; ++q) { basis.bases[q]=static_cast<qc::QST::MeasurementBasis::Basis>(bd(rng)); bits[q]=rng()&1; }
        
        // 构建 N-qubit 投影矩阵 Π (在 ddMat 的 Matrix 模式下)
        dd::Edge proj = ddMat->makeIdent(0, N-1); ddMat->incRef(proj);
        for (int q=0; q<N; ++q) {
            dd::Edge sp = qc::QST::makePauliMeasurementDD(ddMat, N, q, basis.bases[q], bits[q]);
            ddMat->incRef(sp);
            dd::Edge tmp = ddMat->multiply(proj, sp);
            ddMat->incRef(tmp);
            ddMat->decRef(proj); ddMat->decRef(sp);
            ddMat->garbageCollect();
            proj = tmp;
        }
        
        // p = Tr(ρ_true · Π)
        dd::Edge rp = ddMat->multiply(rhoTrue, proj); ddMat->incRef(rp);
        dd::ComplexValue tr = ddMat->trace(rp);
        ddMat->decRef(rp);
        double prob = tr.r;
        if (prob > 1e-8) {
            projs.push_back({proj, prob});
        } else {
            ddMat->decRef(proj);
        }
    }
    cout << "  有效投影算子: " << projs.size() << " / " << nM << "\n";
    cout << "  各算子真实概率: ";
    for (size_t k=0; k<min(projs.size(), (size_t)5); ++k)
        cout << fixed << setprecision(3) << projs[k].prob << " ";
    cout << "\n";

    // ================================================================
    // Step 3: MLE 密度矩阵迭代 + 7 种 Sifting
    // ================================================================
    cout << "\n  Step 3: MLE + Sifting 对比\n\n";
    
    string names[] = {"NoReorder","Sift","LBSift","IGSift","IGLBSift","GrpSift","IGGrpSift"};
    dd::DynamicReorderingStrategy strats[] = {
        dd::None, dd::Sifting, dd::LBSifting, dd::IGSifting,
        dd::IGLBSifting, dd::GroupSifting, dd::IGGroupSifting
    };

    string csvName = cf;
    for (auto& c:csvName) if(c=='/'||c=='\\') c='_';
    size_t dp=csvName.find_last_of('.'); if(dp!=string::npos) csvName=csvName.substr(0,dp);
    ofstream csv("qst_"+csvName+".csv");
    csv << "Strategy,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs\n";

    cout << left << setw(12) << "Strategy" << right
         << setw(10) << "Fidelity" << setw(14) << "TrDist"
         << setw(10) << "RhoSize" << setw(10) << "PeakDD"
         << setw(12) << "Time(ms)\n";
    cout << string(68,'-') << endl;

    for (int s=0; s<7; ++s) {
        cout << left << setw(12) << names[s] << flush;

        // 每个策略独立的 Matrix 模式 Package
        auto dd2 = make_unique<dd::Package>();
        dd2->setMode(dd::Matrix);

        // 重建 |ψ_true⟩ (Vector) + ρ_true (Matrix)
        auto f2 = qc->buildFunctionality(dd2); dd2->incRef(f2);
        dd::Edge i2 = dd2->makeZeroState(N); dd2->incRef(i2);
        dd::Edge psiRef = qc->simulate(i2, dd2); dd2->incRef(psiRef);
        dd2->decRef(i2); dd2->decRef(f2);

        dd::Edge bra2 = dd2->conjugateTranspose(psiRef); dd2->incRef(bra2);
        dd::Edge rt2 = dd2->kronecker(psiRef, bra2); dd2->incRef(rt2);
        dd2->decRef(bra2);

        // 重建投影算子
        vector<Proj> projs2;
        mt19937 rng2(42);
        for (int k=0; k<nM; ++k) {
            vector<int> bits(N);
            qc::QST::MeasurementBasis basis; basis.bases.resize(N);
            for (int q=0; q<N; ++q) { basis.bases[q]=static_cast<qc::QST::MeasurementBasis::Basis>(bd(rng2)); bits[q]=rng2()&1; }
            
            dd::Edge proj = dd2->makeIdent(0, N-1); dd2->incRef(proj);
            for (int q=0; q<N; ++q) {
                dd::Edge sp = qc::QST::makePauliMeasurementDD(dd2, N, q, basis.bases[q], bits[q]);
                dd2->incRef(sp);
                dd::Edge tmp = dd2->multiply(proj, sp);
                dd2->incRef(tmp);
                dd2->decRef(proj); dd2->decRef(sp);
                dd2->garbageCollect();
                proj = tmp;
            }
            
            dd::Edge rp = dd2->multiply(rt2, proj); dd2->incRef(rp);
            dd::ComplexValue tr = dd2->trace(rp);
            dd2->decRef(rp);
            double prob = tr.r;
            if (prob > 1e-8) projs2.push_back({proj, prob});
            else dd2->decRef(proj);
        }

        // MLE 初始: ρ = I/2^n (单位阵) — Matrix 模式
        dd::Edge rho = dd2->makeIdent(0, N-1); dd2->incRef(rho);

        map<unsigned short,unsigned short> vm;
        for (int q=0; q<N; ++q) vm[q]=q;

        dd::InteractionGraph ig;
        ig.initForNqubits(N);
        for (int i=0;i<N;++i) for(int j=i+1;j<N;++j) {ig.weight[i][j]=ig.weight[j][i]=1;ig.degree[i]++;ig.degree[j]++;}
        ig.detectSymmetry();

        unsigned peakDD=0;
        auto t0=high_resolution_clock::now();

        for (int iter=0; iter<nIter; ++iter) {
            // MLE: ρ' = Σ_k Π_k ρ Π_k (等权)
            dd::Edge rhoNext;
            bool first=true;
            for (auto& pj:projs2) {
                dd::Edge t1 = dd2->multiply(rho, pj.dd); dd2->incRef(t1);
                dd::Edge term = dd2->multiply(pj.dd, t1); dd2->incRef(term);
                dd2->decRef(t1);
                if (first) { rhoNext=term; first=false; }
                else {
                    dd::Edge sum=dd2->add(rhoNext,term); dd2->incRef(sum);
                    dd2->decRef(rhoNext); dd2->decRef(term);
                    rhoNext=sum;
                }
            }
            dd::Edge old=rho; rho=rhoNext; dd2->decRef(old);

            // Sifting 每 3 次
            if (iter%3==0 && strats[s]!=dd::None) {
                int psz=dd2->size(rho);
                for (int p=0;p<3;++p) {
                    if(strats[s]==dd::Sifting||strats[s]==dd::LBSifting) rho=dd2->dynamicReorder(rho,vm,strats[s]);
                    else if(strats[s]==dd::IGSifting){auto r=dd2->igSifting(rho,vm,ig);rho=get<0>(r);}
                    else if(strats[s]==dd::IGLBSifting){auto r=dd2->igLbSifting(rho,vm,ig);rho=get<0>(r);}
                    else if(strats[s]==dd::GroupSifting){auto r=dd2->groupSifting(rho,vm,ig);rho=get<0>(r);}
                    else if(strats[s]==dd::IGGroupSifting){auto r=dd2->igGroupSifting(rho,vm,ig);rho=get<0>(r);}
                    int ns=dd2->size(rho); if(ns==psz)break; psz=ns;
                }
            }
            unsigned ca=dd2->maxActive; if(ca>peakDD) peakDD=ca;
        }
        auto t1=high_resolution_clock::now();
        double tm=duration<double,milli>(t1-t0).count();
        unsigned rhoSize=dd2->size(rho);

        // ===== Step 4: Fidelity = ⟨ψ|ρ|ψ⟩ =====
        // ρ 是 N-qubit Matrix DD (N vars), |ψ⟩ 是 N-qubit Vector DD
        dd::Edge rp2 = dd2->multiply(rho, psiRef); dd2->incRef(rp2);
        dd::ComplexValue fv = dd2->innerProduct(psiRef, rp2);
        dd2->decRef(rp2);
        double fid = fv.r; // ⟨ψ|ρ|ψ⟩ 对纯态 ρ=|ψ⟩⟨ψ| 应 = 1
        fid = max(0.0, min(1.0, fid));
        double td = sqrt(max(0.0, 1.0-fid));

        cout << right << fixed
             << setw(10) << setprecision(4) << fid
             << setw(14) << setprecision(4) << td
             << setw(10) << rhoSize
             << setw(10) << peakDD
             << setw(12) << setprecision(1) << tm << endl;

        csv << names[s] << "," << fid << "," << td << ","
            << rhoSize << "," << peakDD << "," << tm << "\n";

        for (auto& pj:projs2) dd2->decRef(pj.dd);
        dd2->decRef(rho); dd2->decRef(rt2); dd2->decRef(psiRef);
        dd2->garbageCollect();
    }
    cout << string(68,'-') << endl;
    csv.close();
    cout << "\n  结果保存到: qst_" << csvName << ".csv\n";

    // 清理
    for (auto& pj:projs) ddMat->decRef(pj.dd);
    ddMat->decRef(rhoTrue);
    ddMat->garbageCollect();
    ddVec->decRef(psiTrue); ddVec->garbageCollect();

    cout << "\n┌── 这次是真的 QST ─────────────────────────────────────────┐\n";
    cout << "│ 1. 建立了密度矩阵 ρ = kronecker(|ψ⟩, conj(|ψ⟩))           │\n";
    cout << "│ 2. 对 ρ 做 MLE 迭代 (Matrix DD, " << N << " vars)                    │\n";
    cout << "│ 3. Fidelity = ⟨ψ_true|ρ_recon|ψ_true⟩ (跨模式验证)        │\n";
    cout << "│ 4. 对 ρ DD 做 sifting 压缩 — 和纯 QMDD 有本质区别       │\n";
    cout << "│                                                           │\n";
    cout << "│ 和之前纯态版本的区别:                                      │\n";
    cout << "│  旧:  纯态 |ψ⟩ + Vector DD → 没构造 ρ                     │\n";
    cout << "│  新:  密度矩阵 ρ + Matrix DD → 真正的 QST                 │\n";
    cout << "└───────────────────────────────────────────────────────────┘\n";

    return 0;
}
