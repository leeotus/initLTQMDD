#!/usr/bin/env python3
"""
compare_qiskit.py — Qiskit Aer 密度矩阵仿真 baseline

与 benchmark_density_matrix (DD方法) 运行相同实验，
输出 CSV 供直接对比：内存占用、仿真时间、物理量正确性。

用法:
  python3 compare_qiskit.py <mode> [options]

mode:
  scale      Dense 内存 vs n（理论值 + 实测仿真时间）
  rounds     多轮噪声下 purity 演化（验证物理正确性）
  expval     <Z> 期望值解析解验证

options:
  --nmin N          最小 qubit 数 (default: 4)
  --nmax N          最大 qubit 数 (default: 10)
  --noise-p P       depolarizing 噪声强度 (default: 0.05)
  --rounds K        rounds 实验轮数 (default: 10)
  --seeds s1,s2     Clifford 种子 (default: 42,123)
  --circuits PATH   外部电路文件或目录（暂不支持，仅合成电路）
  --output FILE     CSV 输出文件 (default: stdout)
"""

import sys
import os

# ---------------------------------------------------------------------------
# 自动安装依赖
# ---------------------------------------------------------------------------
def _ensure_packages():
    required = ["numpy", "qiskit", "qiskit_aer", "psutil"]
    missing = []
    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"[INFO] 安装缺失依赖: {missing}", file=sys.stderr)
        import subprocess
        # 尝试 pip install
        cmd = [sys.executable, "-m", "pip", "install", "--quiet"] + \
              [p.replace("_", "-") for p in missing]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            # pip 本身不存在，尝试 ensurepip
            subprocess.run([sys.executable, "-m", "ensurepip", "--upgrade"],
                           capture_output=True)
            result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("[ERROR] 无法安装依赖，请手动执行:", file=sys.stderr)
            print("  pip install numpy qiskit qiskit-aer psutil", file=sys.stderr)
            sys.exit(1)
        print("[INFO] 依赖安装完成", file=sys.stderr)

_ensure_packages()

# ---------------------------------------------------------------------------
# 正式 import
# ---------------------------------------------------------------------------
import argparse
import time
import tracemalloc
import math

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import QFT as QiskitQFT
from qiskit.quantum_info import (
    DensityMatrix, Kraus, Statevector, SparsePauliOp, random_clifford
)

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------
def _rss_mb():
    """当前进程 RSS（常驻内存）MB"""
    if _HAS_PSUTIL:
        return psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024
    except Exception:
        pass
    return -1.0

def dense_mb(n):
    return 2.0 * (4 ** n) * 8 / (1024 * 1024)

def _depolarizing_kraus(p):
    """构造与 C++ 实现完全一致的 depolarizing Kraus 算子"""
    a = math.sqrt(1.0 - p)
    b = math.sqrt(p / 3.0)
    return Kraus([
        a * np.eye(2, dtype=complex),
        b * np.array([[0, 1], [1, 0]], dtype=complex),
        b * np.array([[0, -1j],[1j, 0]], dtype=complex),
        b * np.array([[1, 0], [0, -1]], dtype=complex),
    ])

def _apply_noise(dm_obj, n, p):
    """对所有 qubit 施加一轮 depolarizing 噪声"""
    channel = _depolarizing_kraus(p)
    for q in range(n):
        dm_obj = dm_obj.evolve(channel, qargs=[q])
    return dm_obj

def _build_qft_dm(n):
    """构建 QFT|0> 的密度矩阵"""
    qc = QuantumCircuit(n)
    qc.compose(QiskitQFT(n, do_swaps=False), inplace=True)
    sv = Statevector.from_instruction(qc)
    return DensityMatrix(sv)

def _build_clifford_dm(n, seed):
    """构建 RandomClifford|0> 的密度矩阵"""
    cliff = random_clifford(n, seed=seed)
    qc = QuantumCircuit(n)
    qc.append(cliff, range(n))
    sv = Statevector.from_instruction(qc)
    return DensityMatrix(sv)

def _purity(dm_obj):
    return np.real(np.trace(dm_obj.data @ dm_obj.data))

def _trace(dm_obj):
    return np.real(np.trace(dm_obj.data))

def _expval_z0(dm_obj, n):
    """计算 Z⊗I⊗...⊗I 期望值（qubit 0 上的 Z）"""
    paulis = "I" * (n - 1) + "Z"  # Qiskit 字符串从高位到低位
    op = SparsePauliOp(paulis)
    return np.real(dm_obj.expectation_value(op))

# ---------------------------------------------------------------------------
# Experiment 1: scale
# CSV: source,type,n,noise_p,dense_nodes,dense_MB,sim_time_ms,peak_mem_MB,purity,trace
# ---------------------------------------------------------------------------
def exp_scale(cfg, out):
    for n in range(cfg.nmin, cfg.nmax + 1):
        builders = [("QFT", lambda n=n: _build_qft_dm(n))]
        for s in cfg.seeds:
            builders.append((f"Clifford_{n}_{s}",
                             lambda n=n, s=s: _build_clifford_dm(n, s)))

        for name, builder in builders:
            typ = "QFT" if "QFT" in name else "Clifford"
            try:
                mem_before = _rss_mb()
                tracemalloc.start()
                t0 = time.perf_counter()

                dm_obj = builder()
                dm_obj = _apply_noise(dm_obj, n, cfg.noise_p)

                sim_ms = (time.perf_counter() - t0) * 1000
                _, peak_bytes = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                mem_after = _rss_mb()

                peak_mb = max(peak_bytes / 1024 / 1024,
                              mem_after - mem_before if mem_after > 0 else 0)
                pur = _purity(dm_obj)
                tr  = _trace(dm_obj)

                out.write(f"{name},{typ},{n},{cfg.noise_p:.3f},"
                          f"{4**n},{dense_mb(n):.2f},"
                          f"{sim_ms:.1f},{peak_mb:.2f},"
                          f"{pur:.6f},{tr:.6f}\n")
                out.flush()
            except Exception as e:
                print(f"skip {name}: {e}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Experiment 2: rounds
# CSV: source,type,n,round,noise_p,dense_nodes,purity,trace,sim_ms_per_round
# ---------------------------------------------------------------------------
def exp_rounds(cfg, out):
    for n in range(cfg.nmin, cfg.nmax + 1):
        builders = [("QFT", lambda n=n: _build_qft_dm(n))]
        for s in cfg.seeds[:1]:  # 只取第一个 seed
            builders.append((f"Clifford_{n}_{s}",
                             lambda n=n, s=s: _build_clifford_dm(n, s)))

        for name, builder in builders:
            typ = "QFT" if "QFT" in name else "Clifford"
            try:
                dm_obj = builder()
                pur = _purity(dm_obj)
                tr  = _trace(dm_obj)
                out.write(f"{name},{typ},{n},0,{cfg.noise_p:.3f},"
                          f"{4**n},{pur:.6f},{tr:.6f},0.0\n")
                out.flush()

                for r in range(1, cfg.rounds + 1):
                    t0 = time.perf_counter()
                    dm_obj = _apply_noise(dm_obj, n, cfg.noise_p)
                    ms = (time.perf_counter() - t0) * 1000
                    pur = _purity(dm_obj)
                    tr  = _trace(dm_obj)
                    out.write(f"{name},{typ},{n},{r},{cfg.noise_p:.3f},"
                              f"{4**n},{pur:.6f},{tr:.6f},{ms:.1f}\n")
                    out.flush()
            except Exception as e:
                print(f"skip {name}: {e}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Experiment 3: expval
# 单 qubit 验证解析解，多 qubit Clifford 展示衰减
# CSV: source,type,n,noise_p,expval_Z0,analytical,purity,trace
# ---------------------------------------------------------------------------
def exp_expval(cfg, out):
    ps = [0.0, 0.01, 0.05, 0.10, 0.20, 0.50]

    # 单 qubit |0> 解析验证
    if cfg.nmin <= 1:
        n = 1
        for p in ps:
            qc = QuantumCircuit(1)  # |0>
            dm_obj = DensityMatrix(Statevector.from_instruction(qc))
            if p > 0:
                dm_obj = _apply_noise(dm_obj, 1, p)
            ev  = _expval_z0(dm_obj, 1)
            ana = 1.0 - 4.0 * p / 3.0
            pur = _purity(dm_obj)
            tr  = _trace(dm_obj)
            out.write(f"|0>,single_qubit,1,{p:.3f},{ev:.6f},{ana:.6f},{pur:.6f},{tr:.6f}\n")
            out.flush()

    # 多 qubit Clifford（首个 seed）
    for n in range(max(cfg.nmin, 2), cfg.nmax + 1):
        s = cfg.seeds[0]
        name = f"Clifford_{n}_{s}"
        for p in ps:
            try:
                dm_obj = _build_clifford_dm(n, s)
                if p > 0:
                    dm_obj = _apply_noise(dm_obj, n, p)
                ev  = _expval_z0(dm_obj, n)
                pur = _purity(dm_obj)
                tr  = _trace(dm_obj)
                out.write(f"{name},Clifford,{n},{p:.3f},{ev:.6f},N/A,{pur:.6f},{tr:.6f}\n")
                out.flush()
            except Exception as e:
                print(f"skip {name} p={p}: {e}", file=sys.stderr)

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Qiskit Aer density matrix baseline")
    parser.add_argument("mode", choices=["scale", "rounds", "expval"])
    parser.add_argument("--nmin",    type=int,   default=4)
    parser.add_argument("--nmax",    type=int,   default=10)
    parser.add_argument("--noise-p", type=float, default=0.05, dest="noise_p")
    parser.add_argument("--rounds",  type=int,   default=10)
    parser.add_argument("--seeds",   type=str,   default="42,123")
    parser.add_argument("--output",  type=str,   default="-")
    args = parser.parse_args()
    args.seeds = [int(s) for s in args.seeds.split(",")]

    out = open(args.output, "w") if args.output != "-" else sys.stdout

    if args.mode == "scale":
        out.write("source,type,n,noise_p,dense_nodes,dense_MB,sim_time_ms,peak_mem_MB,purity,trace\n")
        exp_scale(args, out)
    elif args.mode == "rounds":
        out.write("source,type,n,round,noise_p,dense_nodes,purity,trace,sim_ms_per_round\n")
        exp_rounds(args, out)
    elif args.mode == "expval":
        out.write("source,type,n,noise_p,expval_Z0,analytical,purity,trace\n")
        exp_expval(args, out)

    if args.output != "-":
        out.close()
        print(f"[INFO] 结果已保存: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
