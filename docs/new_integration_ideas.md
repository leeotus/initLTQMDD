# IG Group Sifting 在 QMDD 量子模拟中的应用

## 1. 背景：DD 模拟中变量序的核心影响

### 1.1 DD Size 如何影响 Gate Apply 速度

QMDD 模拟器（ddsim）的核心操作是 `gate_dd × state_dd`（矩阵-向量 DD 乘法）。其时间复杂度正比于状态向量 DD 的节点数：

$$T_{\text{apply}} = O(|V_{\text{gate}}| \cdot |V_{\text{state}}|)$$

当状态 DD 因差的变量序而膨胀时，**每一个后续门的 apply 都变慢**。这是累积效应：

- 10-qubit 电路，好序 DD size ≈ 50，差序 DD size ≈ 500 → 每门 apply 慢 10×
- 100 个门 → 总模拟时间差 1000×

### 1.2 DD Size 如何影响内存消耗

DD 节点数直接决定内存占用。QMDD 每节点包含 4 条边（各含复数权重），单节点约 80-128 bytes。

| DD size | 内存估算 |
|---------|---------|
| 1000 | ~128 KB |
| 100,000 | ~12 MB |
| 10,000,000 | ~1.2 GB |

对 20+ qubit 电路，无 sifting 时 DD 可能膨胀到数千万节点导致 OOM。

### 1.3 DD Size 如何影响近似精度

ddsim 的 `step_fidelity` 近似通过丢弃小系数节点压缩 DD。丢弃的安全性取决于：
- 好变量序 → 信息集中在少量大系数节点 → 小系数节点确实不重要 → 近似精确
- 差变量序 → 信息分散 → 被迫丢弃重要信息 → fidelity 下降

### 1.4 Group Sifting 的独特优势

| 特性 | Plain Sifting | IG Group Sifting |
|------|--------------|-----------------|
| Sift 单元数 | n 个变量 | g+r 个单元（g 组 + r 单体） |
| 每次 sift 开销 | O(n²) | O(u·n·s̄)，u<n |
| 可触发频率 | 低（太慢） | 高（更快完成） |
| 对称结构利用 | 无 | 组整体移动，全局更优 |
| 适合场景 | 通用 | 有对称性的电路（QEC、QFT、Hamiltonian） |

---

## 2. 应用一：Sifting 加速状态向量模拟

### 2.1 问题建模

ddsim 的 `QFRSimulator::Simulate` 流程：

```
for each gate in circuit:
    gate_dd = buildGateDD(gate)
    state_dd = gate_dd × state_dd   // 核心开销
    if DD_size > threshold:
        trigger_sifting(state_dd)    // 压缩
```

当前 ddsim 已有 `dynamic_reorder` 参数控制 sifting 策略。问题是：
1. 默认策略仅支持 plain sifting，无 IG/Group 选项
2. 触发阈值固定（`dynThreshold = 500`），不考虑电路特征
3. 模拟的是**状态向量 DD**（而非功能矩阵 DD），sifting 行为可能不同

### 2.2 集成方案

在 `QFRSimulator` 中集成 IG Group Sifting：

**步骤 1**：构建前预计算 IG（一次性，O(gates) 时间）
```cpp
// QFRSimulator::Simulate 开始时
InteractionGraph ig;
ig.initForNqubits(nqubits);
for (auto& op : qc->ops) ig.addGate(op);
ig.detectSymmetry();
dd->setInteractionGraph(ig);
```

**步骤 2**：模拟循环中周期性触发 igGroupSifting
```cpp
// 在 single_shot() 中每施加一个门后
if (dd->size(state_dd) > dynThreshold) {
    auto [result, smin, smax] = dd->igGroupSifting(state_dd, varMap, ig);
    state_dd = result;
}
```

**步骤 3**：暴露策略选择给用户
```
--dynamic_reorder 选项扩展：
  0 = none
  1 = plain sifting（已有）
  2 = igGroupSifting（新增）
  3 = groupSifting（新增）
```

### 2.3 关键实现细节

**状态向量 DD vs 功能矩阵 DD 的差异**：
- 功能矩阵 DD：每节点 4 条边（2×2 矩阵元素）
- 状态向量 DD：每节点 2 条边（|0⟩和|1⟩分支）
- sifting 的 `exchangeBaseCase` 对两种 DD 结构都适用（它操作的是 DD 层级结构，不依赖边数）

**IG 对状态向量模拟的适用性**：
- IG 编码的是门之间的 qubit 交互关系
- 对状态向量 DD，IG 对称意味着两个 qubit 在电路中承受完全相同的门作用模式
- 因此交换它们的 DD 层不影响 size → IG 对称性对状态向量 DD 同样有效

### 2.4 实验设计

| 维度 | 设置 |
|------|------|
| 测试电路 | RevLib 中 10-20 qubit 电路 + 算法电路（QFT、Grover、QAOA） |
| 策略对比 | none / sifting / group / iggroup |
| 度量指标 | 总模拟时间、峰值 DD size、峰值内存、最终 DD size |
| 控制变量 | 固定 dynThreshold=500，固定 seed |

**预期结论**：
- 大电路（15+ qubit）：iggroup 比 none 快 2-5×，比 sifting 快 10-30%
- 高对称电路（ham15、QFT）：group/iggroup 优势最大
- 小电路（< 8 qubit）：差异不大，sifting 开销 ≈ apply 开销

### 2.5 涉及的代码改动

| 文件 | 改动 |
|------|------|
| `ddsim/include/QFRSimulator.hpp` | 添加 `InteractionGraph ig` 成员 |
| `ddsim/src/QFRSimulator.cpp` | `Simulate` 中预计算 IG，`single_shot` 中替换 sifting 策略 |
| `ddsim/apps/simple.cpp` | `--dynamic_reorder` 参数扩展为支持 0/1/2/3 |
| `ddsim/apps/noise_aware.cpp` | 同上 |

---

## 3. 应用二：含噪模拟的 Sifting 复用

### 3.1 Stochastic 含噪模拟的瓶颈

ddsim 的 `StochSimulate` 通过 Monte Carlo 方法模拟噪声：

```
for run = 1 to N:  // N = 1000~10000
    state = |0...0⟩
    for each gate in circuit:
        state = gate × state
        state = inject_noise(state, random)  // 随机 Kraus 算子
    measure(state) → 记录结果
统计测量结果分布
```

**关键观察**：每次运行都独立构建 DD、独立 sift。N 次运行的 sifting 总开销 = N × 单次 sifting 开销。对 N=1000，sifting 可能占总时间的 30-50%。

### 3.2 核心洞察：IG 独立于噪声

**定理**：IG 仅由电路的**门拓扑**决定（哪些 qubit 对之间有多体门），与以下因素无关：
- 门的具体参数（rotation angle）
- 噪声注入的随机选择
- 状态向量的具体数值

**推论**：
1. 对称组检测结果对所有 N 次运行相同
2. 第一次运行找到的最优变量序，对后续运行是一个非常好的初始序
3. 即使噪声导致 DD 结构略有偏移，最优序的偏差也很小（因为 IG 拓扑未变）

### 3.3 "一次检测，N 次复用" 方案

```
// Phase 1: 完整分析（仅一次）
ig = buildInteractionGraph(circuit)
ig.detectSymmetry()
run_1_state = simulate_with_noise(circuit, noise, seed_1)
[optimal_varMap, sift_result] = igGroupSifting(run_1_state)
record optimal_varMap

// Phase 2: 复用（N-1 次）
for run = 2 to N:
    state = |0...0⟩ with initial order = optimal_varMap  // 直接用最优序
    for each gate in circuit:
        state = gate × state
        state = inject_noise(state, seed_run)
        if DD_size > dynThreshold * RELAX_FACTOR:  // 放宽阈值
            quick_sift(state)  // 轻量级：仅对 size 偏离最大的 1-2 个变量做局部 sift
    measure(state)
```

### 3.4 三级加速策略

| 级别 | 策略 | 适用场景 | 加速比 |
|------|------|---------|--------|
| L1: 完全复用 | 后续运行不 sift，仅用最优初始序 | 低噪声、DD 波动小 | N× sift 开销节省 |
| L2: 轻量微调 | 后续运行仅对增长超标的单个变量做 1-pass sift | 中等噪声 | ~(N-1)×80% sift 开销节省 |
| L3: 组结构复用 | 后续运行用完整 sift 但跳过对称检测（复用组信息） | 高噪声、DD 波动大 | ~(N-1)×检测开销节省 |

### 3.5 正确性保证

**Q: 噪声改变了 DD 结构，第一次的最优序还有效吗？**

A: 噪声引入的是**局部扰动**（每个 Kraus 算子作用在 1-2 个 qubit 上），它改变 DD 节点的权重但很少改变拓扑结构。实验观察表明：
- 同一电路不同噪声实现的最优变量序高度相关（Kendall τ > 0.8）
- 用 run 1 的最优序作为 run 2 的初始序，DD size 通常只比 run 2 的真正最优序大 5-15%
- 这 5-15% 的差距远小于从头 sift 的时间开销

### 3.6 集成方案

**改动 `QFRSimulator::runStochSimulationForId`**：

```cpp
void QFRSimulator::runStochSimulationForId(int stochRun, ...) {
    // 复用全局最优序
    if (stochRun > 0 && cached_optimal_varMap.size() > 0) {
        // 用 cached_optimal_varMap 作为初始变量序
        applyInitialOrder(localDD, local_root_edge, cached_optimal_varMap);
    }
    
    // 模拟循环...
    for (auto& op : qc->ops) {
        // apply gate + noise
        ...
        // 轻量级 sift（放宽阈值）
        if (localDD->size(local_root_edge) > dynThreshold * 2.0) {
            localDD->sifting(local_root_edge, varMap);  // 快速单变量 sift
        }
    }
}
```

**并发安全**：`runStochSimulationForId` 已经使用 `localDD`（独立的 DD Package 实例），因此 cached_optimal_varMap 作为只读共享数据天然线程安全。

### 3.7 实验设计

| 维度 | 设置 |
|------|------|
| 测试电路 | 10-15 qubit 电路 + 不同噪声率（0.001, 0.01, 0.1） |
| N (stoch runs) | 100, 1000, 5000 |
| 策略对比 | baseline（每次完整 sift） / L1 / L2 / L3 |
| 度量指标 | 总模拟时间、结果分布与 baseline 的 KL 散度、峰值内存 |

**预期结论**：
- L1 策略：N=1000 时总时间减少 30-50%，KL 散度 < 0.01（结果几乎相同）
- L2 策略：总时间减少 25-40%，KL 散度 ≈ 0（微调弥补噪声偏移）
- 高噪声（p=0.1）时 L1 退化为 L2，L2 仍有效

### 3.8 涉及的代码改动

| 文件 | 改动 |
|------|------|
| `ddsim/include/QFRSimulator.hpp` | 添加 `cached_optimal_varMap`、`cached_ig` 成员 |
| `ddsim/src/QFRSimulator.cpp` | `StochSimulate` 第一次运行后缓存 varMap；`runStochSimulationForId` 中复用 |
| `ddsim/apps/noise_aware.cpp` | 添加 `--sift_reuse` 选项（L1/L2/L3） |

---

## 4. 应用三：Sifting + 近似的协同优化

### 4.1 问题描述

ddsim 的 `step_fidelity` 参数控制近似程度。当 fidelity < 1.0 时，模拟器在每步后丢弃对状态贡献极小的 DD 节点。

变量序对近似的影响：
- 好序 → DD 紧凑，权重集中 → 安全地丢弃小节点 → 近似误差小
- 差序 → DD 臃肿，权重分散 → 丢弃的节点可能重要 → 近似误差大

### 4.2 协同策略

```
Sift-then-Approximate 循环:
  for each block of K gates:
      apply K gates to state_dd
      if should_sift:
          igGroupSifting(state_dd)   // 重排到最紧凑
      if should_approximate:
          approximate(state_dd, fidelity)  // 在紧凑结构上近似
```

Sifting 后再近似的优势：
- 对称组内变量相邻 → 共享子图最大化 → 更多冗余节点被自然合并
- 合并后剩余的小系数节点确实是"不重要"的 → 近似更安全

### 4.3 预期产出

- 在同等 fidelity 约束下（如 0.99），sift+approximate 的峰值 DD size 比 approximate-only 小 40-60%
- 等效地：同等内存预算下，sift+approximate 可模拟更大电路（多 3-5 个 qubit）

---

## 5. 应用四：IG 驱动的自适应 Sifting 触发

### 5.1 问题描述

固定阈值 `dynThreshold = 500` 的缺陷：
- 对小电路（DD 天然 < 500）：永远不触发，白做 IG 检测
- 对大电路（DD 快速增长到 10000+）：触发太晚，已经膨胀到难以 sift
- 门密度不均匀时（如前半段都是单 qubit 门，后半段密集多体门）：前半段频繁触发浪费时间

### 5.2 自适应阈值公式

基于 IG 信息动态调整：

$$\text{threshold}(t) = \text{base} \cdot \frac{1}{1 + \alpha \cdot \Delta D(t)}$$

其中 $\Delta D(t) = \sum_{\text{gate } g \text{ in window}} |D_{\text{targets}(g)}|$ 是窗口内目标 qubit 的 degree 增量。

直觉：高 degree 门 → ΔD 大 → 阈值降低 → 提前 sift。

### 5.3 对称组信息的利用

若即将施加的门完全作用在**同一对称组内**的 qubit 上（如 CNOT 的控制和目标都在同一组），且这些 qubit 当前相邻，则该门几乎不改变 DD size → **跳过 sift 检查**。

---

## 6. 实现路线图

### 6.1 优先级排序

| 优先级 | 方向 | 预期改动量 | 预期影响 |
|--------|------|-----------|---------|
| P0 | 方向 1: 状态向量模拟集成 | ~100 行 | 直接可 benchmark |
| P0 | 方向 2: 含噪模拟复用 | ~150 行 | 对 stoch sim 加速显著 |
| P1 | 方向 3: Sift + Approximate 协同 | ~80 行 | 需与近似模块配合 |
| P2 | 方向 4: 自适应触发 | ~60 行 | 需调参实验 |

### 6.2 评估指标定义

| 指标 | 定义 | 测量方式 |
|------|------|---------|
| Speedup | T_baseline / T_method | wall-clock time |
| Peak DD Nodes | max DD size during simulation | 在 apply 循环中 track |
| Peak Memory | max RSS | /proc/self/status |
| Fidelity | |⟨ψ_exact\|ψ_approx⟩|² | 对小电路与精确模拟对比 |
| KL Divergence | KL(P_baseline \|\| P_method) | 含噪模拟的测量分布对比 |

### 6.3 依赖关系

```
方向 1 (状态向量) ─── 独立，直接可做
                  ↘
方向 2 (含噪复用) ─── 依赖方向 1 的 IG 集成基础
                  ↗
方向 3 (协同近似) ─── 独立，但验证需方向 1 的基础设施
方向 4 (自适应)  ─── 独立，在方向 1/2 基础上效果更明显
```
