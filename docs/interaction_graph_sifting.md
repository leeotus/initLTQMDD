# Interaction Graph 驱动的 Sifting 变量选择策略

## 1. 动机

标准 QMDD Sifting 算法选择下一个要处理的变量时，仅依据 `active[varMap[j]]`（节点数最多者优先），方向决策仅依据 `pos < n/2`。这一策略忽视了变量（qubit）之间的结构关系——两个频繁交互的 qubit 如果在 DD 层中距离远，会产生大量无法共享的中间节点。

**核心观察**：QMDD 中两个变量层之间的节点共享度取决于对应 qubit 的纠缠/交互程度。Interaction Graph 直接编码了这种结构信息，可以在 sifting 之前就给出哪些变量更"值得优先调整"以及"往哪个方向移动更有利"的先验知识。

---

## 2. Interaction Graph 定义

对量子电路 $C$，构建 **Qubit Interaction Graph** $G = (V, E, w)$：

| 元素 | 定义 |
|------|------|
| $V$ | 节点集 = 电路中所有 qubit $\{q_0, q_1, ..., q_{n-1}\}$ |
| $E$ | 边集 = 所有多 qubit 门涉及的 qubit 对 $\{(q_i, q_j)\}$ |
| $w(q_i, q_j)$ | 边权 = qubit $q_i$ 和 $q_j$ 之间的多 qubit 门数量 |

构建复杂度：$O(|gates| \cdot k^2)$，$k$ = 每门 qubit 数，一次性开销。

---

## 3. IG 增强的两个维度

### 3.1 变量选择：混合评分

标准 sifting 选 `active` 最大的变量。IG 增强后使用混合评分：

$$\text{Score}(q_i) = \alpha \cdot \frac{\text{active}[\text{varMap}[i]]}{\max_j \text{active}[j]} + (1-\alpha) \cdot \frac{C_D(q_i)}{\max_j C_D(q_j)}$$

其中 $C_D(q_i) = \sum_{j \neq i} w(q_i, q_j)$ 为加权度中心性。

**实现中 $\alpha = 0.85$**。实验发现：
- $\alpha$ 过低（如 0.6）时 IG 权重过大，会改变变量处理顺序导致 LB 剪枝误判，结果显著恶化
- $\alpha = 0.85$ 保持 active 节点数为主导，IG 仅作为同等 active 水平下的 tie-breaker

### 3.2 方向决策：邻接引力

标准 sifting 方向：`pos < n/2 → 先下后上`。IG 增强使用引力模型：

$$\text{Gravity}_\text{up}(q_i) = \sum_{j > \text{pos}(q_i)} w(q_i, q_j), \quad \text{Gravity}_\text{down}(q_i) = \sum_{j < \text{pos}(q_i)} w(q_i, q_j)$$

- $\text{Gravity}_\text{down} > \text{Gravity}_\text{up}$ → 先向下 sift（下方有更多交互伙伴）
- 引力相等时退化为标准启发式 `pos < n/2`

---

## 4. 适用范围与实验结论

### 4.1 无剪枝 Sifting（IGSifting）

无剪枝 sifting 中每个变量都完整遍历所有位置，最终找到的最优位置不受处理顺序和方向影响。因此 **IG 对无剪枝 sifting 的最终结果无影响**，仅改变中间路径。

| 电路 | Sifting Size | IGSifting Size | 差异 |
|------|-------------|---------------|------|
| alu4_201 (22q) | 932 | 932 | 无 |

### 4.2 LB Sifting（IGLBSifting）

LB Sifting 有方向剪枝——lower bound 超过当前最优时提前停止某方向。此时方向选择直接决定能否找到更优位置。IG 引力引导先探索更有希望的方向。

当前实验结果（$\alpha = 0.85$）：IG+LB 与纯 LB 结果一致，未产生明显差异。这是因为 $\alpha$ 偏保守，IG 影响被压低。

### 4.3 Linear Sifting（IG + Linear）

IG 同时集成到三种 Linear Sifting 变体中，增强变量选择和方向决策。Linear Sifting 的搜索空间更大（exchange + linear transformation），IG 引导可能在此场景下产生更大差异。

---

## 5. 已实现的策略列表

| 枚举值 | 命令行参数 | 描述 |
|--------|-----------|------|
| `IGSifting` | `ig` | IG 增强的标准 Sifting（无剪枝） |
| `IGLBSifting` | `iglb` | IG 增强的 LB Sifting |
| `igUpperLinearSifting` | `igupperls` | IG 增强的 Upper Linear Sifting |
| `igLowerLinearSifting` | `iglowerls` | IG 增强的 Lower Linear Sifting |
| `igMixLinearSifting` | `igmixls` | IG 增强的 Mix Linear Sifting |

---

## 6. 代码架构

### 6.1 文件结构

```
extern/dd_package/
├── include/
│   ├── InteractionGraph.h      ← IG 数据结构（build/scoring/gravity）
│   └── DDpackage.h             ← 枚举值 + 函数声明
└── src/
    ├── DDigSifting.cpp          ← igSifting + igLbSifting 实现
    ├── DDlinearv2.cpp           ← IG 重载的 linearAndSiftingAux + mixLinearAndSiftingAux
    └── DDreorder.cpp            ← dynamicReorder dispatch

src/
└── QuantumComputation.cpp       ← buildFunctionalityDynamic + parseDynSiftStrategy

apps/
└── main.cpp                     ← 命令行入口
```

### 6.2 InteractionGraph 类

```cpp
struct InteractionGraph {
    int n;                              // qubit 数
    std::vector<std::vector<int>> weight; // weight[i][j] = 门交互次数
    std::vector<int> degree;            // 加权度

    template<typename QC> void build(const QC& qc);  // 从电路构建
    std::vector<short> getSiftOrder() const;          // 按度排序
    std::vector<short> getHybridSiftOrder(...) const; // 混合排序
    bool shouldSiftUpFirst(...) const;                // 引力方向判断
};
```

### 6.3 igSifting / igLbSifting

- `igSifting`：完整复制标准 sifting 逻辑，替换变量选择（hybrid score）和方向决策（gravity）
- `igLbSifting`：完整复制 lbSifting 逻辑，替换变量选择和方向决策，使用 3-arg `exchangeBaseCase` 更新 varMap

### 6.4 IG + Linear Sifting

通过函数重载实现：

```cpp
// 原始版本
Edge linearAndSiftingAux(Edge in, varMap, upOrLow, pruning);

// IG 增强版本（新增 ig 参数）
Edge linearAndSiftingAux(Edge in, varMap, upOrLow, const InteractionGraph& ig, pruning);
```

IG 版本替换变量选择循环和方向判断 `if(pos < n/2)`，其余逻辑（linearSiftingDown/Up、undoMoves、backward）完全复用。

---

## 7. 使用方式

### 7.1 命令行

```bash
# 动态构造 + 最终 IG+LB Sifting
./ltqmdd circuit.real iglb

# 动态构造 + 最终 IG Mix Linear Sifting
./ltqmdd circuit.real igmixls

# 通过环境变量
LTQMDD_DYN_SIFT=igupperls ./ltqmdd circuit.real
```

### 7.2 编程接口

```cpp
qc::QuantumComputation qc("circuit.real");
auto dd = std::make_unique<dd::Package>();

// 方式一：动态构造（构造期间用 sifting 压缩中间态，最终用指定策略优化）
qc::permutationMap map = qc.initialLayout;
auto e = qc.buildFunctionalityDynamic(dd, map, dd::igMixLinearSifting);

// 方式二：先构造后优化
auto e = qc.buildFunctionality(dd);
dd::InteractionGraph ig;
ig.build(qc);
auto [result, minSz, maxSz] = dd->igLbSifting(e, map, ig);

// 方式三：通过 dynamicReorder（IG 为空，退化为标准行为）
auto e = dd->dynamicReorder(e, map, dd::igUpperLinearSifting);
```

### 7.3 Benchmark

```bash
# 对比 Sifting vs LB Sifting vs IG+LB Sifting
./test/run_benchmark_ig_sifting.sh ~/circuits/ results.csv
```

---

## 8. 理论依据

让高权边对应的 qubit 在 DD 中尽量相邻，等价于经典的 **Minimum Linear Arrangement (MLA)** 问题：

$$\min_\pi \sum_{(u,v) \in E} w(u,v) \cdot |\pi(u) - \pi(v)|$$

MLA 是 NP-hard，但 IG 驱动的 sifting 通过贪心策略逼近最优解。

---

## 9. 复杂度分析

| 操作 | 复杂度 | 频率 |
|------|--------|------|
| 构建 IG | $O(\|gates\| \cdot k^2)$ | 一次性 |
| 计算加权度 | $O(n^2)$ | 一次性 |
| 每轮变量选择 | $O(n)$ | 每个变量一次 |
| 方向引力计算 | $O(n)$ | 每个变量一次 |
| 总 IG 额外开销 | $O(n^2)$ | 远小于 sifting 本体 $O(n^2 \cdot \text{exchange\_cost})$ |

---

## 10. 后续优化方向

1. **自适应 $\alpha$**：根据电路特征（门密度、qubit 数）自动调节混合权重
2. **IG + Tight LB**：将 IG 方向引导与更紧的子节点共享度下界结合
3. **归并构造集成**：在 merge-style DD 构造中，用 IG 指导合并前的变量序对齐方向
4. **多轮 IG 更新**：动态构造过程中，根据已处理的门更新 IG 权重（增量式）
