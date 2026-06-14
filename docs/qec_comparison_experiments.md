# QEC 对比实验设计

## 1. 对比工具选择

### 1.1 核心对比工具

| 工具 | 开发者 | 底层方法 | 设计优势 | 设计局限 |
|------|--------|----------|----------|----------|
| **Stim** | Craig Gidney (Google) | Stabilizer Tableau | 万级 qubit Clifford 模拟；Pauli frame noise 传播 | 仅支持 Clifford 门，不支持 T/Toffoli |
| **Qiskit Aer** | IBM | State Vector / Density Matrix | 全门集支持；密度矩阵可精确噪声模拟 | ~30 qubit 状态向量上限 |
| **QuEST** | T. Jones et al. | State Vector (MPI/GPU) | 高性能并行状态向量模拟 | ~30 qubit，无 DD 压缩 |
| **本工具 QEC/QMDD** | 本项目 | QMDD + IG/GroupSifting | DD 压缩 + IG 对称 + 全门集 | 中等规模 (预计 50-100 qubit 对 QEC 电路) |

### 1.2 各工具的安装方式

**Stim** (pip 安装最方便):
```bash
pip install stim
```

**Qiskit Aer**:
```bash
pip install qiskit qiskit-aer
```

**QuEST**:
```bash
git clone https://github.com/QuEST-Kit/QuEST.git
cd QuEST && mkdir build && cd build
cmake .. && make -j$(nproc)
```

**本工具 QEC/QMDD**:
```bash
cd ~/workshop/lbqmdd
cmake -S . -B build_qec -DBUILD_QEC=ON -DGIT_SUBMODULE=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build_qec --target qec_benchmark -j$(nproc)
```

### 1.3 已有论文中的对比基准

以下是已在量子计算模拟/纠错文献中使用的 benchmark 参考：

| 基准来源 | 简介 | GitHub 地址 |
|----------|------|-------------|
| **Stim** (Gidney, Quantum 2021) | 标准 stabilizer 模拟器，自带 benchmark circuits | https://github.com/quantumlib/Stim |
| **QuEST** (Jones et al., SciRep 2019) | 高性能量子模拟器 | https://github.com/QuEST-Kit/QuEST |
| **Qiskit Aer** (IBM, 2019+) | IBM 量子模拟后端 | https://github.com/Qiskit/qiskit-aer |
| **ZX-calculus 等价性检查** (Kissinger & van de Wetering, 2020) | 基于 ZX-diagram 的电路优化 | https://github.com/Quantomatic/pyzx |
| **QEC with Tensor Networks** (Darmawan & Poulin, PRL 2017) | MPS-based Surface Code 模拟 | 无公开代码 (论文方法) |
| **PECOS** (Tannu & Qureshi, ASPLOS 2019) | 噪声感知 qubit 映射 | https://github.com/swamit/tqec |
| **Topological QEC Compiler** (Suchara et al., 2014) | Surface Code 编译 | 算法方法 (概念验证) |

---

## 2. 实验矩阵设计

### 2.1 公共实验设置

```
测试电路集:
  - Steane [[7,1,3]]: encoding + 1-qubit identity (13 physical qubits)
  - Steane [[7,1,3]]: encoding + logical Hadamard (13 physical qubits)
  - 多轮 QEC: 1, 2, 5, 10, 20 轮 syndrome extraction

噪声模型:
  - Depolarization: p ∈ {1e-5, 1e-4, 1e-3, 1e-2}
  - Amplitude Damping + Phase Flip: p ∈ {1e-4, 1e-3}

硬件:
  - 单机: 8-core CPU, 32GB RAM (标准工作站)

评估指标:
  - 模拟时间 (wall-clock)
  - 峰值内存使用
  - 可扩展性 (qubit 数上限)
  - 输出精度 (fidelity / LER 的统计误差)
```

### 2.2 实验对比矩阵

| 实验 | 维度 | Ours (QMDD) | Stim | Qiskit Aer | QuEST |
|------|------|-------------|------|------------|-------|
| A | Clifford 电路正确性 | ✓ | ✓ | ✓ | ✓ |
| B | Clifford QEC 模拟时间 vs 轮数 | ✓ | ✓ | ✓ | ✓ |
| C | 峰值内存 vs 轮数 | ✓ | ✓ | ✓ | ✓ |
| D | Non-Clifford (T 门) QEC 模拟 | ✓ | ✗ | ✓ | ✓ |
| E | DD 压缩率 vs sifting 策略 | ✓ (专属) | ✗ | ✗ | ✗ |
| F | IG 对称组检测 | ✓ (专属) | ✗ | ✗ | ✗ |

---

## 3. 各对比实验的详细设计

### 实验 A：Clifford 电路等效性验证

**目的**：验证 Ours 在 Clifford QEC 电路上的模拟结果与 Stim/Aer 一致。

**电路**：Steane [[7,1,3]], 1 轮 QEC, 1-qubit identity (无噪声)

**步骤**：
1. Ours: 构建电路 → QMDD 模拟 → 1000 shots 测量 → 解码 → 输出逻辑 0/1 分布
2. Stim: 用 stim.Circuit 构建相同电路 → compile_sampler() → 1000 shots → 解析 logical bit
3. Aer: 用 qiskit.QuantumCircuit → AerSimulator → 1000 shots → 解码

**度量**：
- 逻辑错误率 (LER) 应在 0%（无噪声）
- 测量分布应完全一致（Kolmogorov-Smirnov test）

**预期结果**：三者 LER = 0。

### 实验 B：Clifford QEC 模拟时间 vs 轮数

**目的**：对比不同模拟引擎在 Clifford QEC 电路上的性能可扩展性。

**电路**：Steane [[7,1,3]], 1/2/5/10/20/50 轮 QEC, 1-qubit identity (无噪声)

**步骤**：
1. 构建各轮数的物理电路
2. 分别用 Ours (no sifting)、Ours (with sifting)、Stim、Aer 模拟
3. 记录每轮的 wall-clock 时间和峰值内存

**度量**：
- 模拟时间 (秒) vs QEC 轮数
- 峰值内存 (MB) vs QEC 轮数
- 鲁棒性：各工具支持的最大轮数 (内存溢出前)

**预期对比趋势**：

```
Stim:             O(r) 时间增长（stabilizer tableau ~ O(n²) per round）
Aer (statevec):   O(r) 时间增长，但绝对时间高（2^13 状态）
Ours (no sift):   O(r) 时间增长，DD size 随 r 亚线性增长
Ours (sift):      O(r) + sifting overhead，DD size 更紧凑
```

**假设**：Ours 在 qubit 数中等 (13-20) 时比 Aer 快；Stim 在 Clifford 上最快；Ours 的 DD 紧凑性在大量轮次后显示出优势。

### 实验 C：Non-Clifford 门 QEC 模拟能力

**目的**：展示 Ours 在 Non-Clifford QEC 电路上不可替代的价值（Stim 无法处理）。

**电路**：Steane [[7,1,3]], 5 轮 QEC, **逻辑 T 门** (magic state injection)

**步骤**：
1. Ours: QMDD 模拟（天然支持任意酉门）
2. Aer: 状态向量模拟（支持但受限于内存）
3. Stim: 无法处理（Clifford-only）→ 标记为 "N/A"

**度量**：
- 模拟是否可行 (binary)
- 若可行，模拟时间和内存

**预期**：Ours 和 Aer 可以完成，Stim 无法。Ours 预计比 Aer 快（DD 压缩）。

### 实验 D：DD 压缩效果 vs Sifting 策略

**目的**：内部对比——这是我的工具的独特优势。

**电路**：Steane [[7,1,3]], 10 轮 QEC

**对比项**：
- No sifting
- Standard sifting
- IG-guided sifting
- IGGroup sifting

**度量**：
- 最终 DD size (node count)
- 峰值 DD size
- 模拟时间 (含 sifting 开销)

**预期**：IGGroup sifting 实现最大压缩率；Standard sifting 时间最短但压缩较低；No sifting DD size 爆炸。

### 实验 E：噪声阈值验证

**目的**：在 Depolarization 噪声下，验证 Steane Code 的伪阈值 (pseudo-threshold)。

**电路**：Steane [[7,1,3]], 5/10 轮 QEC

**扫描**：p ∈ [1e-5, 1e-3, 1e-2, 5e-2, 1e-1]，Depolarization

**度量**：
- 逻辑错误率 (LER) vs 物理错误率 (p)
- Pseudo-threshold: p* 使得 LER(p*) = p*
- 无编码 baseline: 裸 qubit 的 LER = p

**预期**：Steane pseudo-threshold ≈ 1e-3 到 5e-3 之间（文献值变化大，取决于具体实现细节）。

### 实验 F：IG 对称性分析

**目的**：展示 IG 分析能力对 QEC 电路结构的洞察。

**电路**：Steane [[7,1,3]], Shor [[9,1,3]], [[5,1,3]] (通过继承 QECCode 扩展)

**度量**：
- IG 对称组数量
- 组大小分布
- Group Sifting 压缩率 vs Standard Sifting

**预期**：Steane Code 有 3 个对称组 (2 data + 1 ancilla)，带来 ~3x 压缩增益。

---

## 4. 自动化对比脚本设计

### 4.1 Python 脚本框架

```python
#!/usr/bin/env python3
"""
qec_comparison.py — 自动化 QEC 对比实验脚本
对比: Ours (QMDD) vs Stim vs Qiskit Aer vs QuEST
"""

import subprocess, json, time, os, csv
import numpy as np

# 如果不安装对应的 Python 包，对应测试标记为 SKIP
try: import stim; HAS_STIM = True
except: HAS_STIM = False

try: from qiskit import QuantumCircuit; from qiskit_aer import AerSimulator; HAS_AER = True
except: HAS_AER = False

BENCHMARK_DIR = "qec_comparison_results"
QEC_BENCHMARK_BIN = "./build_qec/qec/qec_benchmark"

# ================== QEC 电路定义 ==================

def build_steane_logical_circuit(logical_gates: list[str]):
    """
    构建 Steane [[7,1,3]] 编码的物理电路
    返回: (standard_op_list, num_physical_qubits)
    """
    # 通过调用 ./qec_benchmark 的导出模式或直接使用 QEC 库
    # 这里用 subprocess 调用定制版本的 qec_benchmark
    pass

# ================== 实验运行函数 ==================

def run_ours(qec_rounds: int, noise_p: float, shots: int = 1000):
    """运行 Ours (QMDD+QEC) 模拟"""
    t0 = time.time()
    result = subprocess.run([QEC_BENCHMARK_BIN, "--rounds", str(qec_rounds),
                              "--noise", str(noise_p), "--shots", str(shots)],
                             capture_output=True, text=True)
    elapsed = time.time() - t0
    # 解析输出得到 LER, DD size, peak memory
    return parse_ours_output(result.stdout, elapsed)

def run_stim(qec_rounds: int, noise_p: float, shots: int = 1000):
    """运行 Stim 模拟"""
    if not HAS_STIM: return {"status": "SKIP (stim not installed)"}
    # Stim 仅支持 Pauli noise (depolarization = Pauli frame noise)
    t0 = time.time()
    c = stim.Circuit()
    # 构建 Steane encoding + syndrome extraction 电路
    c = build_steane_stim_circuit(qec_rounds)
    sampler = c.compile_sampler()
    results = sampler.sample(shots)
    elapsed = time.time() - t0
    ler = compute_logical_error_rate(results)
    return {"time": elapsed, "ler": ler, "shots": shots}

def run_aer(qec_rounds: int, noise_p: float, shots: int = 1000):
    """运行 Qiskit Aer 模拟"""
    if not HAS_AER: return {"status": "SKIP (qiskit-aer not installed)"}
    t0 = time.time()
    qc = build_steane_qiskit_circuit(qec_rounds, noise_p)
    sim = AerSimulator(method='statevector')
    job = sim.run(qc, shots=shots)
    counts = job.result().get_counts()
    elapsed = time.time() - t0
    ler = compute_logical_error_rate_from_counts(counts)
    return {"time": elapsed, "ler": ler, "shots": shots}

# ================== 主实验循环 ==================

def main():
    os.makedirs(BENCHMARK_DIR, exist_ok=True)

    configs = [
        {"rounds": 1,  "noise": 0.0,   "shots": 1000},
        {"rounds": 5,  "noise": 0.001, "shots": 1000},
        {"rounds": 10, "noise": 0.01,  "shots": 500},
    ]

    results = []
    for cfg in configs:
        row = {"rounds": cfg["rounds"], "noise": cfg["noise"]}

        r_ours = run_ours(cfg["rounds"], cfg["noise"], cfg["shots"])
        row.update({f"ours_{k}": v for k, v in r_ours.items()})

        r_stim = run_stim(cfg["rounds"], cfg["noise"], cfg["shots"])
        row.update({f"stim_{k}": v for k, v in r_stim.items()})

        r_aer = run_aer(cfg["rounds"], cfg["noise"], cfg["shots"])
        row.update({f"aer_{k}": v for k, v in r_aer.items()})

        results.append(row)
        print(f"Done: {cfg}")

    # 输出 CSV
    with open(f"{BENCHMARK_DIR}/results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"Results saved to {BENCHMARK_DIR}/results.csv")

if __name__ == "__main__":
    main()
```

### 4.2 运行方式

```bash
# 安装依赖
pip install stim qiskit qiskit-aer numpy

# 启动对比实验
python3 qec_comparison.py
```

### 4.3 预期输出格式

| rounds | noise | ours_time | ours_ler | ours_ddsize | stim_time | stim_ler | aer_time | aer_ler |
|--------|-------|-----------|----------|-------------|-----------|----------|----------|---------|
| 1 | 0.0 | 0.12 | 0.0 | 245 | 0.002 | 0.0 | 0.45 | 0.0 |
| 5 | 0.001 | 0.22 | 0.008 | 1103 | 0.005 | 0.008 | 1.80 | 0.008 |
| 10 | 0.01 | 0.35 | 0.124 | 1980 | 0.010 | 0.124 | 3.60 | 0.124 |

---

## 5. 对比实验的论文呈现

### 5.1 表格：功能对比

| 特性 | Stim | Qiskit Aer | QuEST | Ours (QMDD) |
|------|------|------------|-------|-------------|
| Clifford QEC 模拟 | ✅ | ✅ | ✅ | ✅ |
| Non-Clifford QEC 模拟 | ❌ | ✅ | ✅ | ✅ |
| DD 压缩 | ❌ | ❌ | ❌ | ✅ |
| IG 对称性检测 | ❌ | ❌ | ❌ | ✅ |
| Sifting 优化 | ❌ | ❌ | ❌ | ✅ |
| Density Matrix 模拟 | ❌ | ✅ | ❌ | ❌ (可扩展) |
| GPU/MPI 并行 | ❌ | ✅ | ✅ | ❌ |
| Python API | ✅ | ✅ | ✅ | ⚙️ (C++ + Python 绑定计划) |
| 最大 Clifford qubit 数 | 10000+ | ~30 | ~30 | ~100 (预计, QEC circuit) |

### 5.2 图表：性能 vs 轮数

X 轴: QEC 轮数 (1, 2, 5, 10, 20, 50)
Y 轴: 模拟时间 (log scale)

三条曲线: Stim (bare, fastest), Ours (with sifting), Ours (no sifting), Aer (slowest)

**关键信息**：Stim 最快但在功能上受限；Ours 在中等规模上优于 Aer；Ours 的 DD 压缩优势在大量轮次后趋于明显。

### 5.3 图表：DD 压缩效果

X 轴: Sifting 策略 (None, Standard, IG, IGGroup)
Y 轴: 压缩率 (no-sifting DD size / sifted DD size)

柱状图，展示 IGGroup sifting 的增量收益。

---

## 6. 开源论文项目作为对比参考

以下是与本工具功能相关的**已发表、有开源代码**的论文项目，可作为对比基准：

### 6.1 直接对比（功能重叠）

| 论文 | 代码 | 对比维度 |
|------|------|----------|
| **Stim: a fast stabilizer circuit simulator** (Gidney, Quantum 2021) | https://github.com/quantumlib/Stim | QEC 模拟速度、正确性、可扩展性 |
| **Stochastic Quantum Circuit Simulation Using Decision Diagrams** (Grurl et al., TCAD 2023) | JKQ DDSIM (本项目的上游) | 噪声模拟精度、DD-based 方法扩展对比 |
| **Qiskit Aer: A High-Performance Simulator** | https://github.com/Qiskit/qiskit-aer | 全门集 QEC 模拟能力对比 |
| **Tensor-Network Simulations of the Surface Code** (Darmawan & Poulin, PRL 2017) | 方法对比 (无公开代码) | 方法论对比: DD vs Tensor Network for QEC |

### 6.2 间接对比 / 补充性验证

| 论文 | 代码 | 对比维度 |
|------|------|----------|
| **QuEST** (Jones et al., SciRep 2019) | https://github.com/QuEST-Kit/QuEST | 状态向量 vs DD 方法对比 |
| **PyZX: Large Scale Automated Diagrammatic Reasoning** (Kissinger 2020) | https://github.com/Quantomatic/pyzx | ZX-calculus vs DD 电路优化对比 |
| **Advanced Equivalence Checking for Quantum Circuits** (Burgholzer & Wille, TCAD 2021) | JKQ QFR 框架内 | DD-based 等价性验证方法对比 |
| **Compiling SU(4) Quantum Circuits to IBM QX Architectures** (Zulehner & Wille, ASP-DAC 2019) | JKQ QFR 框架内 | IG-based mapping 方法对比 |
| **Noise-Adaptive Compiler Mappings for Noisy Intermediate-Scale Quantum Computers** (Murali et al., ASPLOS 2019) | https://github.com/Qiskit/qiskit-terra | 噪声感知编译对比 |

### 6.3 QEC 理论参考（不直接对比但需引用）

| 论文 | 说明 |
|------|------|
| Knill, Nature 2005: "Quantum computing with realistically noisy devices" | QEC pseudo-threshold 概念 |
| Gottesman, arXiv:0904.2557 | QEC 基础理论 |
| Raussendorf & Harrington, PRL 2007 | Fault-tolerant quantum computation with high threshold |
| Fowler et al., PRA 2012: "Surface codes: Towards practical large-scale quantum computation" | Surface Code 理论与 threshold |

---

## 7. 实验执行路线

### Phase 1: 基础验证 (Week 1)

- [ ] 安装 Stim, Qiskit Aer（只需 `pip install stim qiskit qiskit-aer`）
- [ ] 在 Cliff 电路上验证 Ours vs Stim vs Aer 的输出一致性
- [ ] 实现 Python 自动化脚本 `qec_comparison.py`

### Phase 2: 性能对比 (Week 1-2)

- [ ] 运行对比实验 A-D
- [ ] 收集所有工具的数据到 CSV
- [ ] 用 matplotlib/seaborn 生成对比图表

### Phase 3: 分析写作 (Week 2-3)

- [ ] 撰写对比分析（与已有方法的具体数值对比）
- [ ] 准备论文的 evaluation section

---

## 参考文献

1. Gidney, C. "Stim: a fast stabilizer circuit simulator." Quantum 5:497, 2021. https://github.com/quantumlib/Stim
2. Grurl, T. et al. "Stochastic Quantum Circuit Simulation Using Decision Diagrams." IEEE TCAD, 2023.
3. Jones, T. et al. "QuEST and High Performance Simulation of Quantum Computers." SciRep 9:10736, 2019. https://github.com/QuEST-Kit/QuEST
4. Darmawan, A. & Poulin, D. "Tensor-Network Simulations of the Surface Code under Realistic Noise." PRL 119:040502, 2017.
5. Kissinger, A. & van de Wetering, J. "PyZX: Large Scale Automated Diagrammatic Reasoning." EPTCS 318, 2020. https://github.com/Quantomatic/pyzx
6. Knill, E. "Quantum computing with realistically noisy devices." Nature 434:39-44, 2005.
7. Gottesman, D. "An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation." arXiv:0904.2557, 2009.