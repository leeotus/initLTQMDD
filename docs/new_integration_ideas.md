# QMDD 应用方向拓展：基于 IG/GroupSifting + ddsim 模拟器的上层场景

## 概述

本项目目前已具备以下核心技术栈：

| 能力 | 描述 |
|------|------|
| **QMDD 紧凑表示** | 以决策图形式高效存储/操作量子酉矩阵和状态向量 |
| **IG 引导重排序** | 基于 qubit 交互图指导 DD 变量序优化 |
| **Group Sifting** | 利用 IG 对称性检测加速变量重排序 |
| **Linear Sifting** | Xor-based 线性变换扩展搜索空间 |
| **Tight LB 剪枝** | 基于子节点共享度的更紧下界加速 sifting |
| **ddsim 模拟** | QMDD-based 状态向量模拟 + 含噪随机模拟 |
| **电路重综合** | DD → 酉矩阵 → Givens 分解 → 门序列 |

基于以上能力，本文档聚焦于**上层应用场景**——即这套技术栈能解决什么实际问题、能产出什么样的研究成果。

---

## 1. 量子编译器的 DD 驱动优化后端

### 应用场景

现有量子编译器（Qiskit/tket/Cirq）的优化 pass 大多基于 gate-level 的规则匹配（如 CNOT cancellation、gate commutation）。这些方法受限于：
- 只能在**局部窗口**内做优化
- 无法感知**全局矩阵结构**
- 编译后的电路仍然存在大量冗余

### 方案：DD-based 全局优化 Pass

将 QMDD 作为编译器的**IR（中间表示）**，在 DD 层面做全局优化：

```
输入电路 → buildFunctionalityDynamic (DD + igGroupSifting) → 紧凑 DD
                                                          ↓
输出电路 ← DD 直接分解 (Shannon/CSD) ← Givens 分析 ← ExtractMatrix
```

**核心价值**：
- Sifting 压缩 DD 等价于**全局重排 qubit 交互顺序**，找到 gate-level 永远无法发现的结构
- IG 编码了电路的**全局交互拓扑**，可指导分解顺序使 CNOT 最少
- 输出电路比输入电路更精简（即使输入已经经过 gate-level 优化）

**产出的研究**：
- 提出 "DD-assisted quantum circuit optimization" 框架
- 在 RevLib/IBM Q 基准集上与 tket/Qiskit O3 对比
- 度量为编译后电路的门数、深度、CNOT 数、在真实硬件的 fidelity

**对标的已有工作**：
- IBM 的 Qiskit transpiler (gate-level)
- Cambridge Quantum 的 tket (peephole optimization)
- Wille 组的 DD-based 等价性检查（但不做编译优化）

---

## 2. 量子纠错码的高效模拟与验证

### 应用场景

量子纠错码（QEC）如 Surface Code、Steane Code、Shor Code 的模拟面临两个挑战：
1. **编码后的 qubit 数爆炸**（1 个逻辑 qubit → 7-17 个物理 qubit）
2. **需要模拟噪声下的纠错过程**（含测量反馈）

### 方案：QEC 电路的 DD 模拟 + 噪声注入

**为什么 DD 适合 QEC**：
- QEC 电路具有**高度规则的结构**（重复的 stabilizer measurement rounds）→ DD 共享率高 → 紧凑
- IG Group Sifting 能自动检测 stabilizer 之间的对称性
- 噪声注入可以直接用 ddsim 的 `noise_effects` 参数（APD）

**具体实现**：
1. 输入：逻辑电路 + QEC 编码方案（如 Steane [[7,1,3]]）
2. 将逻辑电路**自动展开**为编码后的物理电路（logical CNOT → transversal CNOTs）
3. 在 ddsim 中模拟编码后的电路，注入物理级噪声
4. 分析：逻辑错误率 vs 物理错误率 → 验证 error threshold

**产出的研究**：
- 系统地评估不同 QEC 码在各种噪声模型下的逻辑错误率
- 利用 DD 的紧凑性模拟更大码距的 QEC（如 d=5,7 的 Surface Code）
- IG 分析揭示不同 QEC 码的 qubit 交互模式差异

**对标的已有工作**：
- Google 的 Stim（stabilizer simulator，只能做 Clifford 电路）
- Qiskit 的 Aer（全状态向量模拟，受限于 ~30 qubits）
- DD-based 方法可能在 Clifford+non-Clifford 混合电路上找到独特优势

---

## 3. NISQ 设备的噪声感知电路映射（Noise-Aware Qubit Mapping）

### 应用场景

将逻辑量子电路映射到物理设备上（qubit routing）是 NISQ 编译的核心问题。现有方法主要优化**总 CNOT 数**或**电路深度**，很少考虑**设备噪声的不均匀性**。

### 方案：DD + IG 联合驱动的噪声感知映射

**核心洞察**：
- IG 的 `weight[i][j]` 编码了电路对 qubit 对 (i,j) 的"依赖强度"
- 物理设备的噪声模型也可以表示为物理 qubit 之间的"噪声图"（noise graph）
- 最优映射 = 将高依赖的逻辑 qubit 对放在低噪声的物理 qubit 对上

**算法**：
```
1. 构建逻辑电路的 IG: weight[i][j] = qubit pair 交互次数
2. 获取物理设备的噪声特征: noise[a][b] = 物理 qubit pair 的 gate error rate
3. 求解 assignment problem:
   minimize Σ weight[i][j] × noise[π(i), π(j)]
   subject to: π 是 qubit 排列
4. 用 ddsim 模拟映射后的电路，验证 fidelity 改善
```

**关键优势**：DD 模拟可以直接给出精确的**输出 state fidelity**，不需要在真实硬件上运行。

**产出的研究**：
- 提出 "Interaction-Graph-guided noise-aware qubit mapping" 算法
- 对比 IBM Q 设备的默认 mapper（sabre routing）的 fidelity
- 分析不同噪声模型（amplitude damping vs depolarization）对最优映射的影响

**对标的已有工作**：
- Zulehner & Wille (2019) 的 IG-based qubit mapping（但不用噪声信息）
- IBM 的 noise-adaptive compilation
- 本方案结合了 **IG 结构分析 + DD 精确模拟 + 噪声感知** 三者

---

## 4. 量子算法协同设计（Algorithm-Hardware Co-design）

### 应用场景

设计新的量子算法时，需要回答：
- 这个算法在 NISQ 设备上能跑到多大？
- 哪个 qubit layout 最适合这个算法？
- 算法对噪声的鲁棒性如何？

### 方案：基于 DD 的算法快速原型与评估

**工作流**：
```
算法设计 (python)
      ↓
电路生成 (Qiskit/tket)
      ↓
QMDD 构建 + IG 分析 → 评估电路结构的复杂度
      ↓
ddsim 模拟 (无噪声) → 评估理想算法的性能
      ↓
ddsim 模拟 (含噪声) → 评估 NISQ 可行性
      ↓
反馈给算法设计 → 修改算法结构
```

**具体产出**：
- **DD size growth analysis**：绘制算法的 DD size 随 qubit 数/深度的增长曲线，判断可扩展性
- **IG 结构可视化**：用 IG 展示算法的 qubit 交互模式，指导硬件拓扑选择
- **噪声鲁棒性热力图**：扫描噪声参数空间，找出算法的"安全区"

**产出的研究**：
- 对 QAOA/变分算法进行系统的 DD-based 可扩展性分析
- 提出 "DD-complexity" 作为一种新的量子算法复杂度度量
- 对比不同量子算法（Shor/Grover/QFT/QAOA）在相同噪声模型下的表现

---

## 5. 量子电路等价性检查与编译器验证

### 应用场景

编译器的优化 pass 可能引入 bug。如何在**不运行真实量子计算机**的情况下验证编译后的电路与原电路等价？

### 方案：DD-based 等价性检查 + 噪声容忍验证

**已有的基础**：项目已有 `equiv_check.cpp`，通过构建两个电路的 DD 比较。

**扩展方向**：

**a) 大规模等价性检查**
利用 IG 的社区结构将电路分块，分别检查：
- 找到 IG 的连通分量
- 每个分量对应独立的 qubit 子集
- 对每个子集独立检查等价性 → 避免构建巨大 DD

**b) 噪声容忍等价性**
编译后电路在噪声下与原电路"足够接近"就算等价：
- 定义 $\text{FidelityThreshold}$（如 0.99）
- 用 ddsim 的噪声模型模拟两个电路
- 比较输出分布的统计距离

**c) 编译器优化正确性证明**
- 对编译器的每个优化 pass，用 DD 等价性检查验证其正确性
- 自动生成反例（当等价性检查失败时，输出具体的输入状态使两个电路输出不同）

**产出的研究**：
- 系统地验证 Qiskit/tket 的各个优化 pass 的正确性
- 提出 "approximate equivalence checking for NISQ compilers"

---

## 6. 量子机器学习电路的结构分析

### 应用场景

量子机器学习（QML）中，variational circuit（ansatz）的设计是关键。不同的 ansatz 结构导致不同的 expressibility、trainability、entanglement capability。

### 方案：用 IG 和 DD 分析 QML 电路

**IG 分析**：
- IG 的**聚类系数** → 度量电路产生 entanglement 的集中程度
- IG 的**谱半径** → 度量电路的全局连接强度
- IG 的**对称组数量** → 度量 ansatz 的参数冗余度

**DD 分析**：
- DD size 增长曲线 → 度量 ansatz 的 expressibility（太小的 DD = 表达力不足）
- Sifting 后的 DD 压缩率 → 度量 ansatz 的冗余度
- Tight LB 在不同层的值 → 找到 ansatz 的"瓶颈" qubit

**产出的研究**：
- 提出 "DD-based expressibility measure" 替代现有的基于 fidelity 分布的 expressibility
- 对流行的 ansatz（hardware-efficient, QAOA, UCCSD）进行系统的 IG/DD 分析
- 指导设计更高效的 ansatz（DD size 小但 expressibility 高）

---

## 7. 量子密码协议的安全性分析

### 应用场景

量子密钥分发（QKD）、量子数字签名等协议的安全性分析需要模拟**敌手攻击下的量子态演化**。这通常涉及：
- 大维度 Hilbert 空间
- 非标准的攻击电路
- 需要计算 trace distance / fidelity

### 方案：DD-based 安全性分析平台

**核心能力**：
- `dd->partialTrace()` — 追踪掉敌手控制的 qubit，得到诚实方的 reduced state
- `dd->fidelity()` — 比较实际 state 和理想 state 的保真度
- ddsim 噪声模拟 — 评估实际设备的不完美性对安全性的影响

**示例：BB84 QKD 的安全性分析**
```
1. 构建 Alice 制备 + Eve 攻击 + Bob 测量的完整电路
2. 用 ddsim 模拟（可选注入噪声）
3. 计算密钥率 = mutual information(Alice, Bob) - mutual information(Alice, Eve)
4. 扫描 Eve 的不同攻击策略，绘制安全性边界
```

**产出的研究**：
- 用 DD 模拟更大的 QKD 系统（多光子、诱骗态）
- 评估设备无关 QKD 协议在各种噪声模型下的性能
- 提出新的 QKD 安全性证明辅助工具

---

## 8. 量子电路的可解释性分析

### 应用场景

当前的量子电路是"黑盒"。很难回答：
- 哪些 qubit 是"关键路径"？哪些可以去掉？
- 电路的哪一部分贡献了最多的运算？
- 两个看起来不同的电路是否在数学上等价？

### 方案：DD 驱动的电路可视化与解释

**IG 可视化**：用 IG 图展示 qubit 交互模式，一眼看出电路的结构

**DD 层级分析**：每层 (qubit) 的节点数反映了该 qubit 的"信息承载量"

**Sifting 敏感度分析**：
- 对每个 qubit pair，测试交换它们在 DD 中的位置后 DD size 的变化
- 高敏感度对 = 这两个 qubit 的交互对电路结构至关重要
- 低敏感度对 = 这两个 qubit 可以任意交换，电路存在对称性

**噪声敏感度分析**：
- 对每个 gate 注入噪声，用 ddsim 测量输出 fidelity 下降
- 高敏感度 gate = 电路的"薄弱环节"，需要重点保护
- 低敏感度 gate = 可以用近似实现

**产出的研究**：
- 提出 "QMDD-based quantum circuit interpretability" 框架
- 自动生成电路的可解释性报告（关键路径、薄弱环节、冗余部分）
- 与 IBM 的 Qiskit Pulse 可视化对比

---

## 9. 量子-经典混合算法的 DD 加速

### 应用场景

QAOA、VQE 等变分量子算法需要在经典优化器（调整参数）和量子模拟器（评估 cost）之间循环。

### 方案：参数化电路的 DD 复用

**关键挑战**：每轮迭代的参数不同 → 电路不同 → 需要重新构建 DD

**方案 — DD 模板预计算**：
- QAOA 电路是固定结构 + 可变参数：$e^{-i\gamma H_C} e^{-i\beta H_M}$
- 可以**预计算**不依赖参数的部分的 DD（如 mixer 的 DD）
- 运行时只需要将参数化部分与预计算部分相乘
- 利用 `OperationTable` 缓存常用子电路的 DD

**方案 — Sifting 结果复用**：
- QAOA 的参数更新不改变 qubit 交互结构 → IG 不变
- 第一轮的 sifting 结果可以直接复用给后续轮次
- 大幅减少每轮的 sifting 开销

**产出的研究**：
- 提出 "cached DD evaluation for variational quantum algorithms"
- 对比直接状态向量模拟的加速比
- 分析哪些 ansatz 结构最适合 DD 缓存

---

## 10. 量子电路基准测试与性能预测

### 应用场景

给定一个量子电路和一台量子计算机的规格，能否**在运行之前预测**：
- 电路在这台机器上的 fidelity？
- 是否需要 error mitigation？
- 哪个 qubit mapping 最优？

### 方案：DD + IG 驱动的性能预测器

**特征提取**：
- IG 特征：degree 分布、聚类系数、最大交互权重
- DD 特征：DD size、sifting 压缩率、各层 active nodes 分布
- 设备特征：gate error rates、coherence times、connectivity graph

**预测器**：
- 输入：电路 + 设备特征
- 输出：预测的 fidelity、最佳 mapping、推荐 error mitigation 策略
- 模型：简单 ML（随机森林/梯度提升），训练数据来自 ddsim 模拟

**为什么不用真实设备训练**：真实设备太慢、太贵、太不稳定。DD 模拟可以提供**无限量**的精确训练数据。

**产出的研究**：
- 建立 "quantum circuit performance predictor" 模型
- 在 IBM Q 公开数据上验证预测准确度
- 提供 "circuit quality score" 作为 NISQ 电路的 quality metric

---

## 总结：应用方向优先级

| 优先级 | 应用方向 | 所需改动 | 研究价值 | 实现周期 |
|--------|----------|----------|----------|----------|
| ★★★★★ | 量子编译器 DD 优化后端 (§1) | ddsim + circuit_resynth 已有基础 | 高（实用性强） | 2-4 周 |
| ★★★★★ | NISQ 噪声感知映射 (§3) | 新增 qubit mapper + ddsim 验证 | 高（NISQ 刚需） | 4-6 周 |
| ★★★★ | 量子算法协同设计 (§4) | ddsim 已有所有能力 | 中（方法论创新） | 2-3 周 |
| ★★★★ | 等价性检查与编译器验证 (§5) | equiv_check 已有基础 | 高（工具性） | 3-5 周 |
| ★★★ | QEC 高效模拟 (§2) | 新增 QEC 电路生成器 | 高（前沿方向） | 6-10 周 |
| ★★★ | QML 电路结构分析 (§6) | 新增 IG/DD 分析工具 | 中（新兴领域） | 3-6 周 |
| ★★ | 量子密码安全性分析 (§7) | ddsim + partial trace | 中（niche 方向） | 4-8 周 |
| ★★ | 电路可解释性 (§8) | 新增可视化/分析层 | 中（工具性） | 2-4 周 |
| ★★ | 变分算法 DD 加速 (§9) | OperationTable 扩展 | 中（性能优化） | 3-5 周 |
| ★ | 基准测试与性能预测 (§10) | 新增 ML 模型 | 低（验证性） | 4-6 周 |

---

## 关键差异化优势

与直接用 Qiskit/tket/Stim 等工具相比，本 QMDD 框架在这些应用方向上的**独特优势**：

1. **DD 紧凑性** → 能在单机上模拟 30-50 qubit 的电路（远超 state vector 的 30 qubit 上限）
2. **IG 结构分析** → 提供了电路拓扑的数学化表示，可指导智能决策
3. **Sifting 压缩** → 自动发现电路中隐藏的结构冗余
4. **噪声模拟集成** → 在同框架内完成理想模拟和含噪模拟，便于对比
5. **C++ 实现** → 高性能，适合作为编译器后端的底层工具

## 参考文献

[1] Zulehner & Wille, "Compiling SU(4) quantum circuits to IBM QX architectures," ASP-DAC 2019.

[2] Burgholzer & Wille, "Advanced Equivalence Checking for Quantum Circuits," IEEE TCAD 2021.

[3] Hillmich, Zulehner & Wille, "Decision Diagrams for Quantum Computing," Springer 2022.

[4] Grurl et al., "Noise-aware Quantum Circuit Simulation with Decision Diagrams," IEEE TCAD 2023.

[5] Nielsen & Chuang, "Quantum Computation and Quantum Information," Cambridge 2010.

[6] Preskill, "Quantum Computing in the NISQ era and beyond," Quantum 2018.

[7] Cerezo et al., "Variational Quantum Algorithms," Nature Reviews Physics 2021.

[8] Gidney, "Stim: a fast stabilizer circuit simulator," Quantum 2021.

[9] Sivarajah et al., "t|ket⟩: a retargetable compiler for NISQ devices," Quantum Science and Technology 2020.

[10] Cross et al., "Validating quantum computers using randomized model circuits," PRA 2019.