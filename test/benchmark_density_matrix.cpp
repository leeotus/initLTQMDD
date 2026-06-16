/*
 * benchmark_density_matrix.cpp
 *
 * 用法:
 *   benchmark_density_matrix <mode> [options]
 *
 * mode:
 *   scale      DD节点 vs Dense内存，随 n 增长
 *   rounds     多轮噪声 DD size 演化
 *   compress   igGroupSifting/Sifting 对含噪密度矩阵的压缩
 *   expval     <Z> 期望值随噪声衰减（含解析解验证）
 *
 * options:
 *   --nmin N          最小 qubit 数 (default: 4)
 *   --nmax N          最大 qubit 数 (default: 12)
 *   --noise-p P       噪声强度 (default: 0.05)
 *   --rounds K        噪声轮数，用于 rounds 实验 (default: 10)
 *   --seeds S1,S2,..  随机种子列表，用于 Clifford 电路 (default: 42,123)
 *   --circuits DIR    外部电路目录 (.real/.qasm)，与合成电路一起运行
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <memory>
#include <dirent.h>
#include <algorithm>

#include "QuantumComputation.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"
#include "DensityMatrix.hpp"

using namespace std;

// ---------------------------------------------------------------------------
// 参数结构
// ---------------------------------------------------------------------------
struct Config {
    unsigned short nmin    = 4;
    unsigned short nmax    = 12;
    double noiseP          = 0.05;
    int    rounds          = 10;
    vector<unsigned int> seeds = {42, 123};
    string circuitsDir     = "";
};

static Config parseArgs(int argc, char** argv, int startIdx) {
    Config cfg;
    for (int i = startIdx; i < argc; ++i) {
        string a = argv[i];
        if (a == "--nmin"     && i+1<argc) { cfg.nmin    = (unsigned short)stoi(argv[++i]); }
        else if (a == "--nmax"     && i+1<argc) { cfg.nmax    = (unsigned short)stoi(argv[++i]); }
        else if (a == "--noise-p"  && i+1<argc) { cfg.noiseP  = stod(argv[++i]); }
        else if (a == "--rounds"   && i+1<argc) { cfg.rounds  = stoi(argv[++i]); }
        else if (a == "--circuits" && i+1<argc) { cfg.circuitsDir = argv[++i]; }
        else if (a == "--seeds"    && i+1<argc) {
            cfg.seeds.clear();
            string s = argv[++i];
            stringstream ss(s);
            string tok;
            while (getline(ss, tok, ',')) cfg.seeds.push_back((unsigned int)stoul(tok));
        }
    }
    return cfg;
}

// ---------------------------------------------------------------------------
// 电路来源：合成 or 文件
// ---------------------------------------------------------------------------
struct CircuitSource {
    string     name;          // 展示名
    string     type;          // "QFT" / "Clifford" / "file"
    unsigned short nq = 0;
    string     filePath;      // 仅 file 类型
    unsigned int seed = 0;    // 仅 Clifford
};

static vector<CircuitSource> collectSources(const Config& cfg) {
    vector<CircuitSource> sources;

    // 合成电路
    for (unsigned short n = cfg.nmin; n <= cfg.nmax; ++n) {
        sources.push_back({"QFT_" + to_string(n), "QFT", n, "", 0});
        for (unsigned int s : cfg.seeds) {
            sources.push_back({"Clifford_" + to_string(n) + "_" + to_string(s),
                               "Clifford", n, "", s});
        }
    }

    // 外部电路：目录 或 单个文件
    if (!cfg.circuitsDir.empty()) {
        // 辅助函数：读取单文件的 numvars
        auto readNqubits = [](const string& fpath) -> unsigned short {
            ifstream f(fpath);
            string line2;
            while (getline(f, line2)) {
                if (line2.rfind(".numvars", 0) == 0)
                    return (unsigned short)stoi(line2.substr(9));
                if (line2.rfind("qreg", 0) == 0) {
                    auto lb = line2.find('['), rb = line2.find(']');
                    if (lb != string::npos && rb != string::npos)
                        return (unsigned short)stoi(line2.substr(lb+1, rb-lb-1));
                }
            }
            return 0;
        };

        // 判断是目录还是单文件
        DIR* dir = opendir(cfg.circuitsDir.c_str());
        if (dir) {
            // 目录模式：遍历所有 .real/.qasm
            struct dirent* ent;
            vector<string> files;
            while ((ent = readdir(dir)) != nullptr) {
                string fname = ent->d_name;
                if (fname.size() > 5 &&
                    (fname.substr(fname.size()-5) == ".real" ||
                     fname.substr(fname.size()-5) == ".qasm")) {
                    files.push_back(fname);
                }
            }
            closedir(dir);
            sort(files.begin(), files.end());
            for (auto& fname : files) {
                string fpath = cfg.circuitsDir + "/" + fname;
                unsigned short nv = readNqubits(fpath);
                if (nv >= cfg.nmin && nv <= cfg.nmax) {
                    string base = fname.substr(0, fname.rfind('.'));
                    sources.push_back({base, "file", nv, fpath, 0});
                }
            }
        } else {
            // 单文件模式
            ifstream chk(cfg.circuitsDir);
            if (chk.good()) {
                string fpath = cfg.circuitsDir;
                string fname = fpath.substr(fpath.rfind('/') + 1);
                unsigned short nv = readNqubits(fpath);
                if (nv == 0) {
                    cerr << "警告: 无法读取 qubit 数，跳过: " << fpath << "\n";
                } else if (nv < cfg.nmin || nv > cfg.nmax) {
                    cerr << "警告: " << fname << " 有 " << nv
                         << " qubits，不在 [" << cfg.nmin << "," << cfg.nmax
                         << "] 范围内，仍然运行\n";
                    string base = fname.substr(0, fname.rfind('.'));
                    sources.push_back({base, "file", nv, fpath, 0});
                } else {
                    string base = fname.substr(0, fname.rfind('.'));
                    sources.push_back({base, "file", nv, fpath, 0});
                }
            } else {
                cerr << "警告: 无法打开: " << cfg.circuitsDir << "\n";
            }
        }
    }
    return sources;
}

// ---------------------------------------------------------------------------
// 构建含噪 rho (已 incRef，caller 负责 decRef)
// ---------------------------------------------------------------------------
static dd::Edge makeRho(unique_ptr<dd::Package>& pkg,
                        const CircuitSource& src, double p) {
    unsigned short n = src.nq;
    dd::Edge psi;

    if (src.type == "QFT") {
        qc::QFT qft(n);
        pkg->setMode(dd::Vector);
        dd::Edge zero = pkg->makeZeroState(n);
        pkg->incRef(zero);
        psi = qft.simulate(zero, pkg);
        pkg->incRef(psi);
        pkg->decRef(zero);
    } else if (src.type == "Clifford") {
        qc::RandomCliffordCircuit rc(n, n * 2, src.seed);
        qc::permutationMap map;
        for (unsigned short q = 0; q < n; ++q) map[q] = q;
        pkg->setMode(dd::Matrix);
        auto [U, _] = rc.buildFunctionality(pkg, map);
        pkg->incRef(U);
        pkg->setMode(dd::Vector);
        dd::Edge zero = pkg->makeZeroState(n);
        pkg->incRef(zero);
        psi = pkg->multiply(U, zero);
        pkg->incRef(psi);
        pkg->decRef(U); pkg->decRef(zero);
    } else {
        qc::QuantumComputation qc(src.filePath);
        qc::permutationMap map;
        for (unsigned short q = 0; q < n; ++q) map[q] = q;
        pkg->setMode(dd::Matrix);
        auto [U, _] = qc.buildFunctionality(pkg, map);
        pkg->incRef(U);
        pkg->setMode(dd::Vector);
        dd::Edge zero = pkg->makeZeroState(n);
        pkg->incRef(zero);
        psi = pkg->multiply(U, zero);
        pkg->incRef(psi);
        pkg->decRef(U); pkg->decRef(zero);
    }
    pkg->garbageCollect();

    pkg->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(pkg, psi);
    pkg->incRef(rho);
    pkg->decRef(psi);

    if (p > 0) {
        auto kraus = dm::depolarizingKraus(p);
        for (unsigned short q = 0; q < n; ++q) {
            dd::Edge rn = dm::applyKrausChannel(pkg, rho, kraus, q, n);
            pkg->incRef(rn); pkg->decRef(rho); pkg->garbageCollect();
            rho = rn;
        }
    }
    pkg->garbageCollect();
    return rho;
}

static double denseMB(unsigned short n) {
    return 2.0 * (1ULL << (2*n)) * 8.0 / (1024.0*1024.0);
}

// ===========================================================================
// Experiment 1: scale
// CSV: source,type,n,rho_pure_nodes,rho_noisy_nodes,noise_p,dense_MB,dense_nodes,dd_to_dense_ratio
// ===========================================================================
static void expScale(const Config& cfg) {
    for (auto& src : collectSources(cfg)) {
        try {
            unsigned short n = src.nq;
            auto pkg = make_unique<dd::Package>();

            // pure rho
            dd::Edge rhoPure = makeRho(pkg, src, 0.0);
            unsigned int sizePure = pkg->size(rhoPure);

            // noisy rho
            dd::Edge rhoNoisy = makeRho(pkg, src, cfg.noiseP);
            unsigned int sizeNoisy = pkg->size(rhoNoisy);

            unsigned long dNodes = 1ULL << (2*n);
            double dMB = denseMB(n);
            double ratio = (double)dNodes / max(1u, sizeNoisy);

            cout << src.name << "," << src.type << "," << (int)n << ","
                 << sizePure << "," << sizeNoisy << ","
                 << fixed << setprecision(3) << cfg.noiseP << ","
                 << fixed << setprecision(2) << dMB << ","
                 << dNodes << ","
                 << fixed << setprecision(1) << ratio << "\n";

            pkg->decRef(rhoPure);
            pkg->decRef(rhoNoisy);
        } catch (const exception& e) {
            cerr << "skip " << src.name << ": " << e.what() << "\n";
        }
    }
}

// ===========================================================================
// Experiment 2: rounds
// CSV: source,type,n,round,dd_size,purity,trace,dense_nodes
// ===========================================================================
static void expRounds(const Config& cfg) {
    for (auto& src : collectSources(cfg)) {
        // rounds 实验仅对每种 type 的第一个 seed 运行，避免重复
        if (src.type == "Clifford" && src.seed != cfg.seeds[0]) continue;

        try {
            unsigned short n = src.nq;
            auto pkg = make_unique<dd::Package>();
            dd::Edge rho = makeRho(pkg, src, 0.0);  // pure start

            unsigned long dNodes = 1ULL << (2*n);
            auto kraus = dm::depolarizingKraus(cfg.noiseP);

            // round 0
            cout << src.name << "," << src.type << "," << (int)n << ",0,"
                 << pkg->size(rho) << ","
                 << fixed << setprecision(6) << dm::purity(pkg, rho) << ","
                 << fixed << setprecision(6) << pkg->trace(rho).r << ","
                 << dNodes << "\n";

            for (int r = 1; r <= cfg.rounds; ++r) {
                for (unsigned short q = 0; q < n; ++q) {
                    dd::Edge rn = dm::applyKrausChannel(pkg, rho, kraus, q, n);
                    pkg->incRef(rn); pkg->decRef(rho); pkg->garbageCollect();
                    rho = rn;
                }
                cout << src.name << "," << src.type << "," << (int)n << "," << r << ","
                     << pkg->size(rho) << ","
                     << fixed << setprecision(6) << dm::purity(pkg, rho) << ","
                     << fixed << setprecision(6) << pkg->trace(rho).r << ","
                     << dNodes << "\n";
            }
            pkg->decRef(rho);
        } catch (const exception& e) {
            cerr << "skip " << src.name << ": " << e.what() << "\n";
        }
    }
}

// ===========================================================================
// Experiment 3: compress
// igGroupSifting / Sifting / None 对含噪密度矩阵的压缩
// CSV: source,type,n,noise_p,size_noisy,size_ig,size_sift,ig_ratio,sift_ratio,ig_ms,sift_ms
// ===========================================================================
static void expCompress(const Config& cfg) {
    vector<double> ps = {0.0, cfg.noiseP};

    for (auto& src : collectSources(cfg)) {
        for (double p : ps) {
            try {
                unsigned short n = src.nq;
                auto pkg = make_unique<dd::Package>();
                dd::Edge rho = makeRho(pkg, src, p);
                unsigned int sizeNoisy = pkg->size(rho);

                qc::permutationMap vm;
                for (unsigned short q = 0; q < n; ++q) vm[q] = q;

                auto t0 = chrono::high_resolution_clock::now();
                dd::Edge rhoIG = pkg->dynamicReorder(rho, vm, dd::IGGroupSifting);
                double igMs = chrono::duration<double,milli>(
                    chrono::high_resolution_clock::now()-t0).count();
                unsigned int sizeIG = pkg->size(rhoIG);

                for (unsigned short q = 0; q < n; ++q) vm[q] = q;
                auto t1 = chrono::high_resolution_clock::now();
                dd::Edge rhoS = pkg->dynamicReorder(rho, vm, dd::Sifting);
                double siftMs = chrono::duration<double,milli>(
                    chrono::high_resolution_clock::now()-t1).count();
                unsigned int sizeSift = pkg->size(rhoS);

                cout << src.name << "," << src.type << "," << (int)n << ","
                     << fixed << setprecision(3) << p << ","
                     << sizeNoisy << "," << sizeIG << "," << sizeSift << ","
                     << fixed << setprecision(4) << (sizeNoisy>0?(double)sizeIG/sizeNoisy:1.0) << ","
                     << fixed << setprecision(4) << (sizeNoisy>0?(double)sizeSift/sizeNoisy:1.0) << ","
                     << fixed << setprecision(1) << igMs << ","
                     << fixed << setprecision(1) << siftMs << "\n";

                pkg->decRef(rho);
            } catch (const exception& e) {
                cerr << "skip " << src.name << " p=" << p << ": " << e.what() << "\n";
            }
        }
    }
}

// ===========================================================================
// Experiment 4: expval
// <Z_0> 随噪声衰减，单 qubit 含解析解验证
// CSV: source,type,n,noise_p,expval_Z0,analytical,purity,trace
// ===========================================================================
static void expExpval(const Config& cfg) {
    vector<double> ps = {0.0, 0.01, 0.05, 0.10, 0.20, 0.50};

    // 单 qubit 解析验证 (只运行一次)
    if (cfg.nmin <= 1) {
        for (double p : ps) {
            auto pkg = make_unique<dd::Package>();
            pkg->setMode(dd::Vector);
            dd::Edge psi = pkg->makeZeroState(1);
            pkg->incRef(psi);
            pkg->setMode(dd::Matrix);
            dd::Edge rho = dm::densityMatrixFromState(pkg, psi);
            pkg->incRef(rho);
            if (p > 0) {
                auto kraus = dm::depolarizingKraus(p);
                dd::Edge rn = dm::applyKrausChannel(pkg, rho, kraus, 0, 1);
                pkg->incRef(rn); pkg->decRef(rho); pkg->garbageCollect();
                rho = rn;
            }
            vector<dm::PauliTerm> Hz = {{1.0, {3}}};
            double ev = dm::expectationValue(pkg, rho, Hz, 1);
            double ana = 1.0 - 4.0*p/3.0;
            cout << "|0>,single_qubit,1,"
                 << fixed << setprecision(3) << p << ","
                 << fixed << setprecision(6) << ev << ","
                 << fixed << setprecision(6) << ana << ","
                 << fixed << setprecision(6) << dm::purity(pkg, rho) << ","
                 << fixed << setprecision(6) << pkg->trace(rho).r << "\n";
            pkg->decRef(rho);
        }
    }

    // 多 qubit Clifford (首 seed)
    for (auto& src : collectSources(cfg)) {
        if (src.type == "QFT") continue;  // QFT <Z>=0 无对比意义
        if (src.type == "Clifford" && src.seed != cfg.seeds[0]) continue;

        for (double p : ps) {
            try {
                unsigned short n = src.nq;
                auto pkg = make_unique<dd::Package>();
                dd::Edge rho = makeRho(pkg, src, p);

                vector<int> zPaulis(n, 0);
                zPaulis[0] = 3;  // Z on qubit 0
                vector<dm::PauliTerm> Hz = {{1.0, zPaulis}};
                double ev  = dm::expectationValue(pkg, rho, Hz, n);
                double pur = dm::purity(pkg, rho);
                double tr  = pkg->trace(rho).r;

                cout << src.name << "," << src.type << "," << (int)n << ","
                     << fixed << setprecision(3) << p << ","
                     << fixed << setprecision(6) << ev << ","
                     << "N/A" << ","
                     << fixed << setprecision(6) << pur << ","
                     << fixed << setprecision(6) << tr << "\n";
                pkg->decRef(rho);
            } catch (const exception& e) {
                cerr << "skip " << src.name << ": " << e.what() << "\n";
            }
        }
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0]
             << " <scale|rounds|compress|expval>"
             << " [--nmin N] [--nmax N] [--noise-p P] [--rounds K]"
             << " [--seeds s1,s2] [--circuits DIR]\n";
        return 1;
    }
    string mode = argv[1];
    Config cfg  = parseArgs(argc, argv, 2);

    if      (mode == "scale")    expScale(cfg);
    else if (mode == "rounds")   expRounds(cfg);
    else if (mode == "compress") expCompress(cfg);
    else if (mode == "expval")   expExpval(cfg);
    else { cerr << "Unknown mode: " << mode << "\n"; return 1; }
    return 0;
}
