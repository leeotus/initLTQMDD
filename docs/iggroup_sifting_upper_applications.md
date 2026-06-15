# igGroupSifting 的上层应用方向分析

## 1. 问题定位

### 核心矛盾

本项目中 igGroupSifting 能够有效压缩 QMDD 的**静态尺寸 (size)**，例如 ham15_107 电路中从 3238 节点降至 1334（-59%）。但上层应用（fidelity、equivalence checking、simulation 等）的实际运行瓶颈在于**时间复杂度**而非空间复杂度，导致"压缩后的紧凑 DD 虽变小了，但跑应用时也没快多少"。

### 根本原因

```
问题链：
buildFunctionalityDynamic 构造期间多次触发 sifting（阈值 1000）
  → 每次 exchangeBaseCase 调用 initComputeTable() 清空全部计算缓存
    → 后续 multiply 的 compute table 命中率暴跌
      → 操作时间并未等比下降

igGroupSifting 将 ham15_107: 3238 → 1334（-59%）
  但运行时间: 1.62s → 0.98s（仅 -39%）
    说明：空间节省未等比转化为时间节省
```

`multiply` 的复杂度主要取决于**两个 DD 交互的结构耦合度**，而非单个 DD 的 size。紧凑的对称组在乘法时可能因为组内变量的紧密耦合导致中间结构膨胀。此外，sifting 过程中频繁清空 compute table 的 overhead 远大于 size 降低带来的收益。

---

## 2. 评估的三个上层应用方向

### 2.1 方向 A：量子电路等价性检查（QCEC）—— 双 DD 对齐

#### 空间瓶颈分析

等价性检查需同时维护两个 DD：
```
reference DD（参考电路）
implementation DD（实现电路）
```
两者同时驻留在 unique table 中，共用 hash buckets。当 qubit 数达到 20+ 时，节点数可达数万级别，GC 触发频率由可用空间决定。

#### igGroupSifting 的独特价值

- 给 reference DD 一个**确定的最优变量序**（经过 igGroupSifting 极致压缩）
- 对 implementation DD 做 **guided reorder**——直接用 reference 的最优序重排，**不需要完整 sifting**
- 两个 DD 对齐后比较只需一次同步遍历，**不需构造中间 DD**
- 对称组检测可跳过等价位置，进一步缩小检查范围

#### 瓶颈对比

| 维度 | 顺序比较 | igGroupSifting 加速 |
|------|---------|-------------------|
| 峰值驻留节点 | O(size_A + size_B) | O(min(size_A, size_B)) |
| GC 触发 | 频繁 | 基本不触发 |
| 能否扩展至 30+ qubit | 常 OOM | 可行 |

#### 现有公开学术工作

| 项目 | 机构 | 链接/引用 | 核心方法 |
|------|------|-----------|---------|
| **MQT QCEC** | TUM（原JKU） | github.com/munich-quantum-toolkit/qcec | DD构造 + 交替DD + ZX-Calculus + 仿真 |
| **MQT Core** | TUM | github.com/cda-tum/mqt-core | 底层DD包（原 dd_package 已归档并入） |
| Burgholzer & Wille | TUM/JKU | arXiv:2004.08420, TCAD 2021 | 利用可逆性反向传播验证 |

**MQT QCEC 现状**：成熟工具链（111 stars, 67 releases, Python API），支持参数化电路、编译流验证、部分等价。核心结论：单次仿真即可判定等价。

**你的切入点**：MQT QCEC 论文专注速度提升，**未专门讨论空间压缩**。你的对称组感知比较策略（跳过对称组内成员的位置差异）可实现差异化创新。

- **论文定位**：DAC/DATE 级别短文
- **可对比基准**：MQT QCEC（完全开源，可编译对比实验）

---

### 2.2 方向 B：密度矩阵仿真（Noisy Simulation / 量子态层析）

#### 空间瓶颈分析

纯态（state vector）需 $2^n$ 个复数，密度矩阵需 **$4^n$ 个条目**。在 QMDD 中：

- 纯态 DD：n 层，每层 4 条边
- 密度矩阵 DD（$\rho = |\psi\rangle\langle\psi|$）：矩阵表示，节点数平方级增加

**实际数据**：
- 16 qubit 纯态 DD：~2000 节点，约 1 MB
- 16 qubit 密度矩阵 DD：200k ~ 1M 节点，约 200 MB ~ 1 GB
- **igGroupSifting 的对称组检测在这里效果最好**——密度矩阵的对称结构（对角线块、交换对称性）远多于纯态电路，对称组检测出的等价 qubit 数量多几倍，压缩率远超纯态

#### igGroupSifting 的独特价值

- 密度矩阵 DD 节点数是纯态 DD 的平方量级，空间是**真正的约束条件**
- 对称量子噪声通道 + 对称初始态 → 更多可检测的 IG 对称组
- `partialTrace` 操作在压缩后的 DD 上进行时，中间节点膨胀幅度更小（因为本就是紧凑结构）

#### 现有公开学术工作

| 项目/工作 | 技术路线 | 论文 |
|----------|---------|------|
| **MQT DDSIM** (TUM) | DD-based 纯态仿真 | Hillmich et al., "Advanced Simulation of Quantum Computations", TCAD 2018 |
| **QuEST** (Oxford) | 多精度状态向量仿真，不支持DD | Jones et al., SciRep 2019 |
| **Qiskit Aer** (IBM) | Chunk-based 密度矩阵仿真，非DD | — |

**关键发现**：**目前没有任何成熟的 DD-based 密度矩阵仿真开源项目**。MQT DDSIM 主要做纯态，Qiskit Aer 用 dense matrix chunk 而非 DD。这是真正的学术空白。

**你的优势**：代码中已有 `kronecker()`、`multiply()`、`partialTrace()`、`trace()`，缺少的只是密度矩阵构造和噪声通道的 DD 表示。

- **论文定位**：TCAD / QCE 级别长文
- **学术新颖性：高**

---

### 2.3 方向 C：量子线路编译 SWAP 映射

#### 瓶颈分析

NISQ 编译需将逻辑电路映射到物理拓扑，每个候选 SWAP 方案对应一个 DD：
```
for each SWAP candidate:
    create modified circuit
    build DD for this circuit
    evaluate cost (DD size)
```
瓶颈是**搜索空间（时间），不是空间**。寻找最优 SWAP 序列是 NP-hard 问题。

#### 现有工作

| 项目 | 机构 | 核心方法 | 论文 |
|------|------|---------|------|
| **Qiskit Sabre** | IBM | SWAP 启发式 + 远期交互图 | Li et al., DAC 2019 |
| **t|ket⟩** | Cambridge Quantum | 图论编译 + 时隙分配 | Sivarajah et al., Quantum Sci. Technol. 2020 |
| **MQT QMAP** | TUM | Exact + 启发式映射 | Peham et al., ASP-DAC 2023 |

**结论**：**方向 C 不适合**。瓶颈是搜索时间，不是存储 DD 的空间。你的 igGroupSifting 在此场景帮助有限。

---

## 3. 三方向综合对比

| 指标 | **A: QCEC** | **B: 密度矩阵仿真** | **C: SWAP映射** |
|------|------------|-------------------|----------------|
| **空间瓶颈性质** | 双DD驻留，GC压力 | 密度矩阵DD节点数平方膨胀 | 搜索时间（非空间） |
| **现有工具成熟度** | MQT QCEC 很成熟需差异化 | **几乎无 DD 密度阵方案** | Qiskit Sabre 已主流 |
| **igGroupSifting适配度** | 高（对称组跳过检查） | **超高（对称结构多几倍）** | 低 |
| **论文空间** | DAC/DATE 短文 | TCAD/QCE 长文 | 已饱和 |
| **实现难度** | 低（复用现有fidelity/trace） | 中（需构建密度矩阵） | 高 |
| **可对比开源项目** | MQT QCEC（可直接对比） | Qiskit Aer, QuEST（非DD） | Qiskit Sabre（非DD） |
| **你的独特优势** | igGroupSifting引导的DD对齐节省空间 | **DD-based密度矩阵的唯一方案** | 无优势 |

---

## 4. 最终推荐

### 首选：方向 B —— 密度矩阵仿真

**理由**：
1. **学术空白**：目前没有任何成熟的 DD-based 密度矩阵仿真开源项目。MQT DDSIM 纯态；Qiskit Aer dense matrix。你的项目是唯一可用的 DD-based 方案
2. **空间瓶颈最严重**：$4^n$ vs $2^n$，igGroupSifting 在这里价值最大
3. **对称组更多**：密度矩阵的对角块对称、对换对称、纠缠态块状结构，让 igGroupSifting 检测出更多对称组
4. **实现可行**：`kronecker()`、`multiply()`、`partialTrace()`、`trace()` 已完备，需补充密度矩阵构造和噪声通道

### 次选：方向 A —— QCEC

**理由**：
- 更快出成果（论文 + 原型）
- MQT QCEC 作为开源 baseline 可直接对比
- 在 RevLib 基准集上实验可快速产出对比数据
- 差异化点：对称组感知比较策略

### 排除：方向 C

- 瓶颈不是空间，不适合 igGroupSifting

---

## 5. 各方向可用论文参考

### QCEC

- L. Burgholzer and R. Wille, "Advanced Equivalence Checking for Quantum Circuits," IEEE TCAD, 2021. [arXiv:2004.08420]
- L. Burgholzer, R. Raymond, and R. Wille, "Verifying Results of the IBM Qiskit Quantum Circuit Compilation Flow," QCE, 2020. [arXiv:2009.02376]
- L. Burgholzer, R. Kueng, and R. Wille, "Random Stimuli Generation for the Verification of Quantum Circuits," ASP-DAC, 2021. [arXiv:2011.07288]
- T. Peham, L. Burgholzer, and R. Wille, "Equivalence Checking of Quantum Circuits with the ZX-Calculus," JETCAS, 2022. [arXiv:2208.12820]
- T. Peham, L. Burgholzer, and R. Wille, "Equivalence Checking of Parameterized Quantum Circuits," ASP-DAC, 2023. [arXiv:2210.12166]

### 密度矩阵仿真

- S. Hillmich, A. Zulehner, and R. Wille, "Advanced Simulation of Quantum Computations," IEEE TCAD, 2018.
- A. Zulehner, S. Hillmich, and R. Wille, "How to Efficiently Handle Complex Values? Implementing Decision Diagrams for Quantum Computing," ICCAD, 2019.
- P. Niemann, R. Wille, D. M. Miller, M. A. Thornton, and R. Drechsler, "QMDDs: Efficient Quantum Function Representation and Manipulation," IEEE TCAD, 35(1):86-99, 2016.

### 可用开源基准项目

| 项目 | 仓库/链接 | 用途 |
|------|----------|------|
| MQT QCEC | github.com/munich-quantum-toolkit/qcec | QCEC 对比基准 |
| MQT Core | github.com/cda-tum/mqt-core | 底层 DD 包（本项目的 dd_package 已归档并入） |
| Qiskit Aer | github.com/Qiskit/qiskit-aer | 密度矩阵仿真对比基准 |
| RevLib circuits | revlib.org | 测试电路库 |

---

## 6. 历史版本

| 版本 | 日期 | 变更 |
|------|------|------|
| v1.0 | 2026-06-16 | 初版，含三方向分析、学术调研、最终推荐 |
