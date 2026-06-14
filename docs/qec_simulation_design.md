# 量子纠错码（QEC）的高效 DD 模拟与验证：详细设计

## 目录

1. [问题定义](#1-问题定义)
2. [QEC 基础回顾](#2-qec-基础回顾)
3. [为什么 DD 特别适合 QEC](#3-为什么-dd-特别适合-qec)
4. [系统架构设计](#4-系统架构设计)
5. [核心模块设计](#5-核心模块设计)
6. [关键算法：QEC 电路自动展开](#6-关键算法qec-电路自动展开)
7. [关键算法：DD 构建中的 IG 对称性利用](#7-关键算法dd-构建中的-ig-对称性利用)
8. [噪声模型与纠错流程模拟](#8-噪声模型与纠错流程模拟)
9. [输出指标与分析](#9-输出指标与分析)
10. [实验设计](#10-实验设计)
11. [实现路线图](#11-实现路线图)
12. [与已有工具的对比定位](#12-与已有工具的对比定位)

---

## 1. 问题定义

### 1.1 核心问题

给定一个**逻辑量子电路** $C_L$（作用于 $k$ 个逻辑 qubit）和一个**QEC 编码方案** $\mathcal{E}$（如 Steane [[7,1,3]]、Shor [[9,1,3]]、Surface [[d²,1,d]]），我们希望回答：

1. **编码后的物理电路** $C_P$ 在噪声信道 $\mathcal{N}$ 下的输出保真度是多少？
2. **逻辑错误率** $p_L$ 与**物理错误率** $p$ 的关系是什么？（error threshold 验证）
3. 不同的 QEC 编码在相同的噪声模型下，哪个表现更好？
4. 纠错 cycle 数对逻辑错误率的影响如何？

### 1.2 当前的瓶颈

| 工具 | 瓶颈 |
|------|------|
| Stim (Gidney) | 仅支持 Clifford 电路（stabilizer formalism），无法处理 T 门等非 Clifford 门 |
| Qiskit Aer | 全状态向量模拟，~30 qubit 上限；QEC 展开后通常 >30 qubit |
| QuEST/Qulacs | 同上，且无 QEC 专用优化 |
| 解析方法 | 仅适用于简单噪声模型和特定编码；无法处理任意逻辑电路 |

### 1.3 DD-based 方法的独特优势

- **DD 的图结构压缩** 天然适配 QEC 电路的**高度规则性**
- **IG + Group Sifting** 能自动检测 stabilizer 之间的**结构对称性**
- **增量 IG** 随纠错 cycle 重复逐步积累，保持对称检测准确性
- 可以同时处理 Clifford+non-Clifford 门，不限于 stabilizer formalism

---

## 2. QEC 基础回顾

### 2.1 Stabilizer 码的结构

一个 $[[n, k, d]]$ stabilizer 码由以下要素定义：

- **$n$ 个物理 qubit**：编码后的 qubit
- **$k$ 个逻辑 qubit**：信息 qubit
- **$d$ 为码距**：可纠正 $\lfloor (d-1)/2 \rfloor$ 个任意错误
- **$n-k$ 个 stabilizer generators**：由 Pauli 算子的张量积构成，用于检测错误
- **逻辑算子** $\bar{X}_i, \bar{Z}_i$：等效于对第 $i$ 个逻辑 qubit 做 X 或 Z 操作

### 2.2 纠错流程（Syndrome Extraction）

每一轮纠错 cycle 包含：

```
1. 噪声作用于物理 qubit（每个 qubit 独立的 error channel）
2. 对每个 stabilizer generator S_j：
   a. 准备一个辅助 qubit |0⟩
   b. 以 S_j 的张量积结构为控制模式，对辅助 qubit 做受控门
   c. 测量辅助 qubit → 得到 syndrome bit
3. 根据 syndrome bit 的组合（syndrome pattern）推断错误位置
4. 施加纠正操作（Pauli X 和/或 Z）
```

### 2.3 目标 QEC 编码

| 编码 | 参数 | 特点 |
|------|------|------|
| **Steane Code** | [[7,1,3]] | 基于经典 Hamming [7,4,3] 码；stabilizer 具有高度对称性；所有 Clifford 门 transversal |
| **Shor Code** | [[9,1,3]] | 级联编码；9 qubit 的重复结构；stabilizer 分组清晰 |
| **Five-Qubit Code** | [[5,1,3]] | 最小可能的纠错码；stabilizer 结构紧凑 |
| **Surface Code (rotated)** | [[d²,1,d]] | 规模化候选；需要大量 qubit；stabilizer 是局部 plaquette/star 算子 |

**初始实现重点**：从 Steane [[7,1,3]] 和 Shor [[9,1,3]] 开始。

---

## 3. 为什么 DD 特别适合 QEC

### 3.1 规则结构 → 高压缩率

QEC 电路的 stabilizer extraction 具有**高度重复的结构**：

```
对每个 stabilizer generator S_j:
  重复相同的模板：
    |+⟩_ancilla --H--●--●--H-- Measure
                     |  |
    data_j1   -----X--|-- ...
    data_j2   --------X-- ...
```

这种重复意味着 QMDD 的唯一表（Unique Table）会有极高的命中率。在 Gate-Level Simulator 中，每个门独立计算；在 DD 中，重复结构自然共享。

**定量示例**：$r$ 轮 syndrome extraction，每轮 $m = n-k$ 个 stabilizer，每个 stabilizer 平均作用于 $w$ 个数据 qubit，总共 $\approx r \cdot m \cdot w$ 个 CNOT 门。但 DD 节点数可能只是 $O(r \cdot n)$ 而非指数级。

### 3.2 IG 对称性检测 ⇔ Stabilizer 代数对称性

**关键观察**：stabilizer 码的 IG 具有以下对称性质：

- **所有数据 qubit** 参与的 stabilizer 完全相同 → **IG 对称组**
- **所有 syndrome qubit** 参与的 stabilizer 测量完全相同 → **IG 对称组**
- Steane Code 中，7 个数据 qubit 分为**2 组对称群**（X-stabilizer 和 Z-stabilizer 对称模式不同）

```
Steane [[7,1,3]] IG 对称组示例：

X-stabilizer group: {q0, q1, q3}  → 3 个 qubit 参与 X 检测模式完全相同
                     {q2, q4, q5, q6} → 4 个 qubit 参与 X 检测模式完全相同

Z-stabilizer group: {q0, q2, q4}  → 3 个 qubit 参与 Z 检测模式完全相同
                     {q1, q3, q5, q6} → 4 个 qubit 参与 Z 检测模式完全相同
```

**Group Sifting 的直接收益**：每组只需 sift 一个代表，其余成员跟随放置。Steane Code 的 7 个数据 qubit 只需要 sift 4 个代表（而非 7 个）。

### 3.3 增量 IG 与多轮纠错的天然契合

```python
for round in range(num_rounds):
    for gate in syndrome_extraction_round:  # 固定模板
        qmdd.apply(gate)
        incremental_ig.addGate(gate)  # IG 逐步累积

    # 每轮结束后，IG 反映该轮的稳定结构
    # Symmetry detection 在每轮都是一致的
```

增量 IG 在多轮重复后自动收敛到稳定的对称组，且不会因为后续 round 的相同门而改变检测结果（幂等性）。

### 3.4 噪声模拟的 DD 优势

噪声算子是**局部的、单 qubit 或两 qubit 的**。在 DD 中：
- 单 qubit 噪声（APD/Pauli）的 DD 表示只有 $O(1)$ 个节点
- 在 sifting 最优的变量序下，噪声算子作用于顶层 qubit，DD 乘法 cost 最小
- 多个 noise operator 的叠加（Kraus sum）可以通过 DD 的 add 操作紧凑表示

---

## 4. 系统架构设计

### 4.1 整体架构

```
┌────────────────────────────────────────────────────────────┐
│                        用户输入                              │
│  · QEC 编码 (Steane/Shor/Surface/... )                      │
│  · 逻辑电路 (QASM/Real)                                     │
│  · 噪声模型 (APD + depolarization + phase flip 参数)         │
│  · 纠错轮数 r                                                │
│  · 模拟模式 (精确/含噪/近似)                                 │
└──────────────────────────┬─────────────────────────────────┘
                           │
                           ▼
┌────────────────────────────────────────────────────────────┐
│                   QEC 电路生成器 (新增)                       │
│  · encode_logical_circuit(logical_circ, code_scheme)        │
│  · generate_syndrome_extraction_round(code_scheme)          │
│  · generate_error_correction(code_scheme, syndrome)          │
└──────────────────────────┬─────────────────────────────────┘
                           │
                           ▼ (物理电路)
┌────────────────────────────────────────────────────────────┐
│              QMDD 构建 + 模拟引擎 (复用+扩展)                 │
│  · buildFunctionalityDynamic (DD + IG 增量)                 │
│  · Group Sifting (IG 对称性检测)                            │
│  · ddsim noise simulation (已支持 APD/X/Y/Z)                │
│  · 新增: syndrome-based conditional correction              │
└──────────────────────────┬─────────────────────────────────┘
                           │
                           ▼
┌────────────────────────────────────────────────────────────┐
│                    输出与分析模块                             │
│  · Logical Error Rate vs Physical Error Rate                │
│  · DD size growth curve vs rounds                           │
│  · IG symmetry analysis report                              │
│  · Final state fidelity heatmap                             │
└────────────────────────────────────────────────────────────┘
```

### 4.2 新增文件

```
qec/
├── include/
│   ├── QECCode.hpp          // QEC 编码的抽象基类
│   ├── SteaneCode.hpp       // Steane [[7,1,3]] 实现
│   ├── ShorCode.hpp         // Shor [[9,1,3]] 实现
│   ├── FiveQubitCode.hpp    // [[5,1,3]] 实现
│   ├── SurfaceCode.hpp      // Surface Code (参数化)
│   └── QECSimulator.hpp     // QEC 模拟器 (继承或组合 QFRSimulator)
├── src/
│   ├── SteaneCode.cpp
│   ├── ShorCode.cpp
│   ├── FiveQubitCode.cpp
│   ├── SurfaceCode.cpp
│   └── QECSimulator.cpp
├── apps/
│   ├── qec_benchmark.cpp     // QEC 基准测试主程序
│   └── CMakeLists.txt
└── test/
    ├── test_qec.cpp
    └── CMakeLists.txt
```

---

## 5. 核心模块设计

### 5.1 QECCode 基类

```cpp
class QECCode {
public:
    // 编码参数
    virtual int n_physical() const = 0;   // 物理 qubit 数
    virtual int k_logical() const = 0;     // 逻辑 qubit 数
    virtual int distance() const = 0;      // 码距

    // Stabilizer 信息
    virtual std::vector<std::vector<PauliOp>> stabilizer_generators() const = 0;
    // S_j[i] = PauliOp (I/X/Y/Z) 作用在第 i 个物理 qubit 上

    // 逻辑算子
    virtual std::vector<PauliOp> logical_X() const = 0;
    virtual std::vector<PauliOp> logical_Z() const = 0;

    // 核心方法：将逻辑门展开为物理门序列
    virtual std::vector<Gate> encodeSingleQubitGate(
        const Gate& logical_gate, int logical_qubit
    ) const = 0;

    virtual std::vector<Gate> encodeTwoQubitGate(
        const Gate& logical_gate, int logical_q1, int logical_q2
    ) const = 0;

    // 生成一轮 syndrome extraction 的电路
    virtual QuantumCircuit generateSyndromeExtraction() const = 0;

    // 根据 syndrome 生成纠错操作
    virtual std::vector<Gate> generateCorrection(
        const std::vector<int>& syndrome
    ) const = 0;
};
```

### 5.2 SteaneCode 实现要点

Steane [[7,1,3]] 的 stabilizer generators:

```
X-stabilizers (测量 X 错误):
  S1 = X I X I X I X   (qubits 0,2,4,6)
  S2 = I X X I I X X   (qubits 1,2,5,6)
  S3 = I I I X X X X   (qubits 3,4,5,6)

Z-stabilizers (测量 Z 错误):
  S4 = Z I Z I Z I Z   (同 S1 但用 Z)
  S5 = I Z Z I I Z Z   (同 S2 但用 Z)
  S6 = I I I Z Z Z Z   (同 S3 但用 Z)
```

**编码逻辑 CNOT 门**：Steane Code 支持 **transversal CNOT**——即对配对物理 qubit 逐对做 CNOT：
```
logical_CNOT(qL_control, qL_target)
→
for i in 0..6:
    physical_CNOT(phys_control[i], phys_target[i])
```

这使得编码后的 CNOT 门**保持电路深度为 1**，不会产生额外开销。

**编码逻辑 T 门**：不是 transversal，需要用**魔法态注入**（magic state injection）+ **gadget** 实现。初始版本可以不支持 T 门，先用 Clifford 电路验证框架。

### 5.3 Syndrome Extraction 电路生成

```python
def generate_syndrome_extraction_round(code):
    """
    生成一轮 syndrome extraction 电路

    对每个 stabilizer generator S_j:
      1. 分配一个辅助 qubit，初始化为 |0⟩
      2. 如果 S_j 包含 X:
         - Hadamard 辅助 qubit (使其变成 |+⟩)
         - 对 S_j 中每个 X/Y 位置，从辅助 qubit 到数据 qubit 做 CNOT
         - Hadamard 辅助 qubit
      3. 如果 S_j 包含 Z:
         - 对 S_j 中每个 Z/Y 位置，从数据 qubit 到辅助 qubit 做 CNOT
      4. 测量辅助 qubit → syndrome bit

    返回: (电路, 辅助 qubit 到 stabilizer index 的映射)
    """
```

**总 qubit 数**：$7 + 6 = 13$ qubit（Steane: 7 data + 6 ancilla），远在 DD 可处理范围内。

### 5.4 QECSimulator

```cpp
class QECSimulator {
    std::unique_ptr<QECCode> code;
    std::unique_ptr<qc::QuantumComputation> physical_circuit;
    std::unique_ptr<QFRSimulator> ddsim;

    // 核心模拟流程
    struct QECResult {
        double logical_error_rate;
        double output_fidelity;
        std::vector<double> round_fidelities;
        unsigned long max_dd_size;
        unsigned long total_dd_nodes;
        std::vector<int> syndrome_history;
    };

    QECResult simulate(
        const qc::QuantumComputation& logical_circuit,
        const NoiseModel& noise,
        int num_rounds,
        int num_shots
    );
};
```

---

## 6. 关键算法：QEC 电路自动展开

### 6.1 主算法

```
Algorithm: EncodeLogicalCircuit

Input:  logical_circuit C_L (k logical qubits, g gates)
        QEC code E ([[n, k, d]])
        num_syndrome_rounds r

Output: physical_circuit C_P (n + ancilla physical qubits)

Procedure:
  C_P ← empty circuit with n + m physical qubits
       (m = n-k syndrome ancillae)

  // Step 1: Encode initial state
  C_P.append( encode_logical_|0⟩^k(E) )

  // Step 2: Interleave logical gates with QEC rounds
  for each gate G in C_L:
      C_P.append( encode_gate(G, E) )

      // Insert syndrome extraction after every t gates
      if gate_count % t == 0:
          C_P.append( generate_syndrome_extraction(E) )
          C_P.append( conditional_error_correction(E) )

  // Step 3: Final syndrome extraction
  C_P.append( generate_syndrome_extraction(E) )
  C_P.append( conditional_error_correction(E) )

  // Step 4: Decode (measure logical qubits)
  C_P.append( decode_logical_state(E) )

  return C_P
```

### 6.2 初始态编码

Steane Code 的 $|0\rangle_L$ 准备：

```
|0⟩_L = 1/√8 (|0000000⟩ + |1010101⟩ + |0110011⟩ + |1100110⟩
             + |0001111⟩ + |1011010⟩ + |0111100⟩ + |1101001⟩)

用电路实现:
  1. 初始化 7 个 qubit 为 |0⟩
  2. 对 qubit 0,1,3 做 Hadamard
  3. 对 (0,2), (0,4), (1,2), (1,5), (3,4), (3,5), (2,6), (4,6), (5,6) 做 CNOT
  → 即用 classical Hamming [7,4,3] 编码的保护方式
```

在 QMDD 中，这个编码电路会产生特定的 DD 结构。关键观察：**编码后的 $|0\rangle_L$ 的 DD size 往往小于任意态的 DD size**，因为编码后的态具有特殊对称性。

### 6.3 逻辑门编码（以 Steane 为例）

| 逻辑门 | 物理实现 | 电路深度 | CNOT 数 |
|--------|----------|----------|---------|
| X | transversal X on data qubits | 1 | 0 |
| Z | transversal Z on data qubits | 1 | 0 |
| H | transversal H on each data qubit | 1 | 0 |
| S | transversal S† on each data qubit (注意 dagger) | 1 | 0 |
| CNOT | transversal CNOT pairs | 1 | 7 |
| T | 需要 magic state + teleportation gadget | ~50 | ~30 |

**关键**：前 5 种门都是 transversal → 电路深度恒为 1 → DD 构建开销与逻辑 qubit 数无关。

---

## 7. 关键算法：DD 构建中的 IG 对称性利用

### 7.1 Multi-Round IG 演化的幂等性

**性质**：在多轮 QEC 中，每一轮 syndrome extraction 的门序列**完全相同**。因此：

- 第 1 轮后的 IG = $G^{(1)}$
- 第 2 轮添加相同门后，$G^{(2)} = 2 \times G^{(1)}$（所有权重翻倍）
- IG 对称组的划分**不变**（仅权重缩放）

这意味着 **Group Sifting 可以在第 1 轮后检测对称组，后续所有轮次复用**。

### 7.2 QEC-Aware Sifting 策略

```
Algorithm: QEC-Aware Sifting Trigger

// 不按固定阈值触发，而是：
// - 每轮 syndrome extraction 开始前做一次 sifting
// - 每轮 syndrome extraction 内部不做 sifting

for round in 1..r:
    // Sifting at round boundary
    dd = igGroupSifting(dd, varMap, ig)
    dd = garbageCollect(dd)

    // Build round without sifting interruptions
    for gate in syndrome_extraction_round:
        dd = multiply(gate_dd, dd)
        incrementalIG.addGate(gate)
```

**为什么这样做**：
- 一轮内的 gate 序列固定且已知良好排序
- 频繁 sifting 在重复电路中是浪费（每轮都一样）
- 轮间 sifting 足以响应噪声导致的微小变化

### 7.3 对称组预计算

对于已知的 QEC 码，对称组可以**离线预计算**而非在线检测：

```cpp
// Steane [[7,1,3]] 的预计算对称组
const QECSymmetry STEANE_SYMMETRY = {
    .n_data_qubits = 7,
    .groups = {
        { /* X-type symmetry */
            {0, 1, 3},    // 参与 X-stabilizer 的对称组 α
            {2, 4, 5, 6}  // 参与 X-stabilizer 的对称组 β
        },
        { /* Z-type symmetry */
            {0, 2, 4},    // 参与 Z-stabilizer 的对称组 γ
            {1, 3, 5, 6}  // 参与 Z-stabilizer 的对称组 δ
        }
    },
    .ancilla_symmetry = { /* 6 个 ancilla 全对称 */ }
};
```

用预计算组替代运行时 `detectSymmetry()`，开销从 $O(n^2)$ 降为 $O(1)$。

---

## 8. 噪声模型与纠错流程模拟

### 8.1 噪声模型

ddsim 已支持的噪声类型：

| 噪声 | 代码标记 | 物理含义 |
|------|----------|----------|
| Amplitude Damping | `A` | $T_1$ 弛豫（能量衰减） |
| Depolarization | `D` | 各向同性误差 |
| Phase Flip | `P` | $T_2$ 退相干（相位误差） |
| Bit Flip (X) | 可组合接入 | 随机 X 错误 |
| Phase+Bit (Y) | 可组合接入 | 随机 Y 错误 |

### 8.2 噪声注入策略

```cpp
// 在每个物理 gate 后以概率 p 注入噪声
for (auto& gate : physical_circuit) {
    dd_apply(gate);

    // 对 gate 涉及的每个 qubit，以概率 p 注入单 qubit 噪声
    for (auto& qubit : gate.getInvolvedQubits()) {
        if (random() < noise_prob) {
            applyNoiseOperator(qubit, noise_type);
        }
    }
}

// 额外：idle qubit 也遭受噪声（一个纠错 cycle 内未参与门的 qubit）
for (auto& qubit : idle_qubits) {
    if (random() < idle_noise_prob) {
        applyNoiseOperator(qubit, noise_type);
    }
}
```

### 8.3 Syndrome-Driven 条件纠错

**挑战**：真实 QEC 需要根据 syndrome measurement 结果条件性地施加纠正操作。但标准 DD 不支持条件分支（measurement feedback）。

**方案 1**：Stochastic unraveling（已有 ddsim 支持）
```
for shot in 1..num_shots:
    full_circuit_with_random_noise = generate_noisy_circuit()
    result = ddsim.simulate(full_circuit_with_random_noise)
    syndrome = extract_syndrome(result)
    correction = lookup_correction(syndrome)
    result_corrected = apply_correction(result, correction)
```

**方案 2**：DD 层面的 direct noise channel（避免 stochastic）
```
ρ_encoded = U_encode |ψ⟩⟨ψ| U_encode†
for round in 1..r:
    ρ = N(ρ)                // noise channel (Kraus sum as DD)
    ρ = S(ρ)                // syndrome extraction (unitary)
    ρ = C_syndrome(ρ)       // conditional correction (weighted sum of possibilities)
output = decode(ρ)
```

**推荐**：先实现方案 1（利用已有 ddsim 能力），方案 2 作为后续优化。

### 8.4 Stochastic 模拟中的复用优化

```cpp
// Build the noiseless DD ONCE
auto dd_perfect = buildPerfectDD(physical_circuit_without_noise);

for (int shot = 0; shot < num_shots; shot++) {
    // Clone the DD (or use separate Package instance per shot)
    auto dd_shot = clone(dd_perfect); // 浅拷贝 + 引用计数

    // Apply random noise operators for this shot
    for (auto& noise_op : generate_random_noise_ops(noise_prob)) {
        dd_shot = dd_shot.apply(noise_op);
    }

    // Measure
    results[shot] = dd_shot.measure();
}
```

**收益**：DD 构建只做一次（最昂贵的操作），每 shot 只做噪声应用 + 测量。

---

## 9. 输出指标与分析

### 9.1 核心指标

| 指标 | 定义 | 计算方法 |
|------|------|----------|
| **Logical Error Rate (LER)** | $p_L = \mathbb{P}[ \text{解码后逻辑态不正确} ]$ | 模拟 N shots，统计解码错误的比例 |
| **Output Fidelity** | $F = \langle \psi_{\text{ideal}} | \rho_{\text{noisy}} | \psi_{\text{ideal}} \rangle$ | ``dd->fidelity(ideal_edge, noisy_edge)`` |
| **DD Size** | DD 节点总数 | ``dd->size(root_edge)`` |
| **DD Compression Ratio** | (sifting 前 size) / (sifting 后 size) | Group Sifting 前后的 size 比 |
| **Pseudo-Threshold** | $p^*$ 使 $p_L = p$ | 用二分搜索在 LER 曲线上找到 |
| **IG Symmetry Count** | 检测到的对称组数 | InteractionGraph.symmetricGroups.size() |

### 9.2 分析报告示例

```
=== QEC Simulation Report ===
Code:            Steane [[7,1,3]]
Logical Circuit: QFT_3 (3 logical qubits, 18 gates)
Physical Qubits: 7 data + 6 ancilla = 13 total
Noise Model:     Depolarization, p = 0.001
QEC Rounds:      10
Shots:           10000

--- Results ---
Logical Error Rate:   0.0234 ± 0.0030
Output Fidelity:      0.9872
Pseudo-Threshold:     0.0085 (estimated)
Max DD Size:          8471 nodes
DD Compression:       3.2x (2650 nodes after Group Sifting)
IG Symmetry Groups:   3 (2 data groups + 1 ancilla group)
Total Simulation Time: 4.2s

--- Round-by-Round Fidelity ---
Round 0: 1.0000
Round 1: 0.9987  ↓
Round 2: 0.9973  ↓
Round 3: 0.9959  ↓
...
Round 10: 0.9872
```

---

## 10. 实验设计

### 10.1 实验 1：Error Threshold 验证

**目标**：验证 Steane Code 的 pseudo-threshold，与理论值对比。

**设置**：
- QEC: Steane [[7,1,3]]、Shor [[9,1,3]]、[[5,1,3]]
- 逻辑电路：单 qubit identity（最简情况，纯测 QEC 效率）
- 噪声：depolarization，扫描 $p \in [10^{-5}, 10^{-1}]$
- 纠错轮数：1, 2, 5, 10

**预期输出**：$p_L(p)$ 曲线，pseudo-threshold $p^*$ 在 $p_L = p$ 处。

**对比基准**：
- Steane Code 理论 pseudo-threshold ~ 0.001-0.01（与具体实现有关）
- 无编码的裸 qubit 错误率曲线

### 10.2 实验 2：DD 可扩展性分析

**目标**：展示 DD 方法在 QEC 模拟中的可扩展性。

**设置**：
- 逐步增加 qubit 数：Steane (13 qubit) → Shor (19 qubit) → concatenated (≥ 49 qubit)
- 记录 DD size vs native state vector size
- 记录模拟时间 vs qubit 数

**预期**：DD size 与 qubit 数呈**亚指数**关系（对于 QEC 电路），而 state vector 是指数关系。

### 10.3 实验 3：逻辑电路保真度

**目标**：评估 QEC 保护下运行实际算法的增益。

**设置**：
- 逻辑电路：QFT_2、GHZ_2、CNOT gate + T gate
- QEC：Steane [[7,1,3]]
- 噪声：APD + Depolarization
- 对比：编码模拟 vs 直接在噪声物理 qubit 上运行逻辑电路

**预期**：当 $p < p^*$ 时，编码版本输出保真度 > 未编码版本。

### 10.4 实验 4：IG 对称性分析

**目标**：量化 QEC 电路的结构对称性，验证 Group Sifting 的收益。

**设置**：
- 分析 IG 的对称组数量、组大小分布
- 对比 Group Sifting vs 标准 Sifting 的 DD size 和速度
- 不同编码的对称性模式对比

**预期**：QEC 码的对称组数 $\approx$ 不同交互模式的 qubit 种类数（通常 2-4 组）。

---

## 11. 实现路线图

### Phase 1: 基础框架（1-2 周）

- [ ] 实现 `QECCode` 抽象基类和 `SteaneCode`
- [ ] 实现 `encode_logical_circuit()` 展开
- [ ] 实现 `generate_syndrome_extraction()` 
- [ ] 编写 `qec_benchmark.cpp` 主程序骨架
- [ ] 单元测试：验证编码后的电路在无噪声下与原电路等价

### Phase 2: 噪声模拟集成（1-2 周）

- [ ] 将 ddsim 噪声模型应用于 QEC 电路
- [ ] 实现 noise-injection-per-gate
- [ ] 实现 idle-qubit noise
- [ ] 实现 stochastic unraveling 模拟
- [ ] 输出 LER、Fidelity 指标

### Phase 3: DD 优化（1 周）

- [ ] 实现 QEC 对称组预计算（离线检测）
- [ ] Multi-round IG 幂等性利用（仅在轮间 sifting）
- [ ] 实现 per-shot DD 复用（noiseless DD 构建一次）
- [ ] 性能 profiling 和优化

### Phase 4: 扩展与发布（1-2 周）

- [ ] 实现 `ShorCode` 和 `FiveQubitCode`
- [ ] 编写实验脚本（自动扫描参数、生成图表）
- [ ] 对比 Stim/Qiskit Aer 结果
- [ ] 撰写文档和 paper outline

---

## 12. 与已有工具的对比定位

| 维度 | Stim | Qiskit Aer | Qulacs | **本 QMDD+QEC** |
|------|------|------------|--------|------------------|
| **电路类型** | Clifford only | All | All | **All** |
| **qubit 上限** | 10000+ | ~30 | ~30 | **~50-100 (QEC 电路)** |
| **噪声模型** | Pauli frame | APD/Depol | Kraus | **APD/Depol/Phase (已有)** |
| **QEC 专用优化** | Stabilizer 快速模拟 | 无 | 无 | **IG 对称组 + Group Sifting** |
| **DD 压缩** | 无 | 无 | 无 | **有 (核心能力)** |
| **等价性检查** | 无 | 无 | 无 | **有 (dd->fidelity)** |
| **开箱即用** | ✅ | ✅ | ❌ (需编程) | ⚙️ (正在构建) |

**核心定位**：本工具填补了 **Clifford+non-Clifford QEC 电路的高效模拟** 这一空白。Stim 在纯 Clifford 上极快但无法处理 T 门；Aer 功能全面但受限于 state vector；本方案借助 DD 的紧凑性在中间地带（50-100 qubit、含非 Clifford 门）找到优势。

---

## 13. 已有工作的文献调研

> 本节回答：**是否有已有的 QMDD（或 Decision Diagram）作用于 QEC 的工作？**

### 13.1 核心结论

**在 Decision Diagram 用于量子纠错码模拟这一具体问题上，目前没有发现已发表的直接工作。** 具体来说：

- **QMDD 文献**（Niemann et al. 2016, Zulehner & Wille TCAD 2018, Hillmich et al. 2022 Springer）聚焦于 DD 的表示能力、构建效率、变量重排序和通用量子电路模拟，**未涉及 QEC 专用模拟**
- **DD-based 噪声模拟**（Grurl et al. TCAD 2023，即 ddsim 的前身）支持通用噪声模型但未做 QEC 特化
- **量子纠错模拟工具**（Stim, Qiskit Aer, QuEST）不使用 DD 作为底层表示

因此，**将 QMDD（含 IG/Group Sifting 等重排序技术）应用于 QEC 模拟是一个新颖的交叉方向**。

### 13.2 直接相关的 DD/量子工作

#### 13.2.1 QMDD 的奠基工作

**P. Niemann, R. Wille, D. M. Miller, M. A. Thornton, R. Drechsler. "QMDDs: Efficient Quantum Function Representation and Manipulation." IEEE TCAD, 35(1):86–99, 2016.**

- QMDD 的定义、唯一性证明、基础操作（add, multiply, Kronecker）
- 实验在 RevLib 可逆电路集上验证（非 QEC 电路）
- **未涉及量子纠错码**

#### 13.2.2 DD Package 与通用模拟

**A. Zulehner, S. Hillmich, R. Wille. "How to Efficiently Handle Complex Values? Implementing Decision Diagrams for Quantum Computing." ICCAD 2019.**

- 描述 JKQ DD Package（本项目的上游）的设计与实现
- 重点在复数缓存（Complex Table）、计算表优化
- **未涉及 QEC 场景**

#### 13.2.3 DD-based 噪声感知模拟（核心前驱）

**T. Grurl, R. Kueng, J. Fuß, R. Wille. "Stochastic Quantum Circuit Simulation Using Decision Diagrams." IEEE TCAD, 42(1):307–321, 2023.**
**（前身为 DATE 2021 会议论文 "Considering Decoherence Errors in the Simulation of Quantum Circuits Using Decision Diagrams"）**

- **这是与本项目最直接相关的工作**——本项目的 ddsim 即基于该工作
- 在 DD 模拟框架中引入噪声算子（Amplitude Damping, Phase Flip, Depolarization）
- 使用 Stochastic unraveling 方法处理噪声通道
- 实验在标准量子算法上验证（Grover, Shor, QFT），**不含 QEC 电路**
- 噪声注入方法**与 QEC 的 syndrome extraction 流程正交**——前者是通用的噪声模拟，后者需要 syndrome 反馈的条件纠错

#### 13.2.4 DD-based 等价性检查

**L. Burgholzer, R. Wille. "Advanced Equivalence Checking for Quantum Circuits." IEEE TCAD, 40(9):1810–1824, 2021.**

- 用 DD 验证两个量子电路的等价性
- 提出了截断和近似等价检查
- **可用于验证编码后电路与原逻辑电路的等价性**（无噪声下），但不涉及噪声模拟或纠错

#### 13.2.5 DD-based 编译与映射

**A. Zulehner, R. Wille. "Compiling SU(4) Quantum Circuits to IBM QX Architectures." ASP-DAC 2019, pp. 185–190.**

- 使用 Interaction Graph (IG) 指导 qubit 映射（逻辑 qubit → 物理 qubit）
- IG 编码了 qubit 间的交互频率，用于求解最小 SWAP 插入
- **IG 分析与 QEC 的分析需求高度重叠**——QEC 的 stabilizer 结构会在 IG 上产生清晰的对称模式

#### 13.2.6 DD 在量子计算中的综述

**S. Hillmich, A. Zulehner, R. Wille. "Decision Diagrams for Quantum Computing." In: Design Automation of Quantum Computers, Springer, pp. 1–26, 2022.**

- 综述了 DD 在量子计算中的各类应用：模拟、等价检查、编译、验证
- **未提及 QEC 作为单独的应用场景**

### 13.3 相邻领域的相关工作

#### 13.3.1 Tensor Network 用于 QEC 模拟

**A. Darmawan, D. Poulin. "Tensor-Network Simulations of the Surface Code under Realistic Noise." PRL 119:040502, 2017.**

- 使用 Tensor Network (MPS/PEPS) 模拟 Surface Code 的纠错过程
- Tensor Network 和 DD 在数学上存在对应关系（两者都是 tensor 的低秩分解）
- **DD 相对于 Tensor Network 的优势**：DD 能精确表示非低秩结构（如 T 门的 magic state），而 TN 在非低秩场景下 bond dimension 会爆炸

#### 13.3.2 Stim — Stabilizer 快速模拟

**C. Gidney. "Stim: a fast stabilizer circuit simulator." Quantum 5:497, 2021.**

- 使用 stabilizer tableau 做 Clifford 电路的快速模拟
- 支持 noise channel（通过 Pauli frame propagation）和 syndrome extraction
- **局限**：不支持非 Clifford 门（T, Toffoli, 任意 angle rotation）
- **本项目的优势**：QMDD 天然支持任意酉门，对于含 T 门的 QEC 电路（如容错实现）有不可替代的价值

#### 13.3.3 Qiskit Aer — 全状态向量模拟

**Qiskit Development Team. "Qiskit Aer: A High-Performance Simulator for Quantum Circuits." 2021.**
**（基于 "Qiskit: An Open-source Framework for Quantum Computing", doi:10.5281/zenodo.2562110, 2019）**

- 支持 state vector, density matrix, stabilizer 等多种模拟方法
- **密度矩阵模拟可用于 QEC**，但矩阵大小随 qubit 数指数增长
- DD 的紧凑表示在处理具有大量结构的 QEC 电路时有潜力超越密度矩阵方法

#### 13.3.4 QuEST / Qulacs — 高性能量子模拟器

**T. Jones et al. "QuEST and High Performance Simulation of Quantum Computers." Scientific Reports 9:10736, 2019.**
**Y. Suzuki et al. "Qulacs: a fast and versatile quantum circuit simulator for research purpose." Quantum 5:559, 2021.**

- 纯状态向量模拟，针对多核/GPU 优化
- 不支持 DD 压缩，状态向量大小始终为 $2^n$
- **当 n~30 时不再是优势**

#### 13.3.5 BDD 用于可逆逻辑的错误检测（概念相关）

**R. Wille, R. Drechsler. "Effect of BDD Optimization on the Synthesis of Reversible and Quantum Logic." IET Electronics Letters, 2012.**
**S.-A. Wang, C.-Y. Lu, I.-M. Tsai, S.-Y. Kuo. "An XQDD-Based Verification Method for Quantum Circuits." IEICE Trans. Fundamentals, 2008.**

- XQDD (X-decomposition Quantum Decision Diagram) 用于量子电路验证
- 利用 BDD-based 方法做可逆电路的错误检测（如 missing gate faults）
- **与 QEC 不同**：这是测试电路制造缺陷，非运行时错误纠正

### 13.4 本项目的创新点定位

综合上述文献分析，**"QMDD + IG/Group Sifting + ddsim 噪声模拟 → QEC 高效模拟"这个方向有明确的创新空间**：

| 维度 | 已有工作 | 本项目计划 |
|------|---------|-----------|
| **DD 用于 QEC** | **无** | **首次探索** |
| DD + 噪声模拟 | Grurl et al. 2023 (通用) | 扩展到 QEC 专用的 syndrome-based 反馈 |
| IG 对称分析 | Zulehner & Wille 2019 (mapping) | 用于检测 stabilizer 对称性 |
| Group Sifting 加速 | 本项目已有工具 | 应用于 QEC 电路的特殊对称结构 |
| Clifford+non-Clifford QEC 模拟 | Stim 不支持非 Clifford | QMDD 天然支持 |

**最好的发表对标**：
- 若在 **TCAD/DAC/DATE** 发：强调 DD-based QEC simulator 的工具性和效率
- 若在 **Quantum/PRX Quantum** 发：强调 "novel simulation methodology enabling the study of QEC codes with non-Clifford gates under realistic noise"
- 若在 **ICCAD/ASP-DAC** 发：强调 DD 优化的技术细节（对称性利用、Group Sifting 加速）

### 13.5 推荐的 Paper 标题与抽象框架

**建议标题**（几个选项）：
1. "Harnessing Decision Diagram Symmetries for Efficient Simulation of Quantum Error-Correcting Codes"
2. "QMDD-based Simulation of Quantum Error Correction: Exploiting Stabilizer Structure via Interaction Graphs"
3. "Efficient QEC Circuit Simulation Using Interaction-Graph-Guided Decision Diagram Reordering"

**核心 claim**：
> We demonstrate that QMDD with IG-based symmetry detection (Group Sifting) can simulate QEC-protected quantum circuits with non-Clifford gates at scales beyond state-vector simulators, by exploiting the structural repetition and stabilizer symmetries inherent to QEC codes — a synergy that has not been explored in prior DD or QEC literature.

---

## 参考文献（完整列表）

[1] D. Gottesman, "An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation," arXiv:0904.2557, 2009.

[2] A. Steane, "Error correcting codes in quantum theory," PRL 77(5):793–797, 1996.

[3] P. Shor, "Scheme for reducing decoherence in quantum computer memory," PRA 52(4):R2493–R2496, 1995.

[4] C. Gidney, "Stim: a fast stabilizer circuit simulator," Quantum 5:497, 2021.

[5] A. Cross et al., "Validating quantum computers using randomized model circuits," PRA 100:032328, 2019.

[6] E. Knill, "Quantum computing with realistically noisy devices," Nature 434:39–44, 2005.

[7] S. Hillmich, A. Zulehner, R. Wille, "Decision Diagrams for Quantum Computing," in Design Automation of Quantum Computers, Springer, pp. 1–26, 2022.

[8] P. Niemann, R. Wille, D. M. Miller, M. A. Thornton, R. Drechsler, "QMDDs: Efficient Quantum Function Representation and Manipulation," IEEE TCAD, 35(1):86–99, 2016.

[9] A. Zulehner, S. Hillmich, R. Wille, "How to Efficiently Handle Complex Values? Implementing Decision Diagrams for Quantum Computing," ICCAD 2019.

[10] T. Grurl, R. Kueng, J. Fuß, R. Wille, "Stochastic Quantum Circuit Simulation Using Decision Diagrams," IEEE TCAD, 42(1):307–321, 2023. (前身: DATE 2021)

[11] A. Zulehner, R. Wille, "Compiling SU(4) Quantum Circuits to IBM QX Architectures," ASP-DAC 2019, pp. 185–190.

[12] L. Burgholzer, R. Wille, "Advanced Equivalence Checking for Quantum Circuits," IEEE TCAD, 40(9):1810–1824, 2021.

[13] A. Darmawan, D. Poulin, "Tensor-Network Simulations of the Surface Code under Realistic Noise," PRL 119:040502, 2017.

[14] T. Jones et al., "QuEST and High Performance Simulation of Quantum Computers," Scientific Reports 9:10736, 2019.

[15] Y. Suzuki et al., "Qulacs: a fast and versatile quantum circuit simulator for research purpose," Quantum 5:559, 2021.

[16] R. Wille, R. Drechsler, "Effect of BDD Optimization on the Synthesis of Reversible and Quantum Logic," IET Electronics Letters, 2012.

[17] S.-A. Wang, C.-Y. Lu, I.-M. Tsai, S.-Y. Kuo, "An XQDD-Based Verification Method for Quantum Circuits," IEICE Trans. Fundamentals, 2008.

[18] A. Zulehner, R. Wille, "Advanced Simulation of Quantum Computations," IEEE TCAD, 2018.
