# Interaction Graph 驱动的 Sifting 变量选择策略

## 1. 动机

标准 QMDD Sifting 算法选择下一个要处理的变量时，仅依据 `active[varMap[j]]`（节点数最多者优先）。这一策略忽视了变量（qubit）之间的结构关系——两个频繁交互的 qubit 如果在 DD 层中距离远，会产生大量无法共享的中间节点。

**核心观察**：QMDD 中两个变量层之间的节点共享度取决于对应 qubit 的纠缠/交互程度。Interaction Graph 直接编码了这种结构信息，可以在 sifting 之前就给出哪些变量更"值得优先调整"以及"往哪个方向移动更有利"的先验知识。

---

## 2. Interaction Graph 定义

对量子电路 $C$，构建 **Qubit Interaction Graph** $G = (V, E, w)$：

| 元素 | 定义 |
|------|------|
| $V$ | 节点集 = 电路中所有 qubit $\{q_0, q_1, ..., q_{n-1}\}$ |
| $E$ | 边集 = 所有多 qubit 门涉及的 qubit 对 $\{(q_i, q_j)\}$ |
| $w(q_i, q_j)$ | 边权 = qubit $q_i$ 和 $q_j$ 之间的多 qubit 门数量 |

---

## 3. 中心性指标

### 3.1 加权度中心性 (Weighted Degree Centrality)

$$C_D(q_i) = \sum_{j \neq i} w(q_i, q_j)$$

含义：该 qubit 与其他 qubit 的**总交互次数**。值高 → 位置敏感度高。

### 3.2 邻接引力 (Adjacency Gravity)

$$\text{Gravity}_\text{up}(q_i) = \sum_{j > \text{pos}(q_i)} w(q_i, q_j), \quad \text{Gravity}_\text{down}(q_i) = \sum_{j < \text{pos}(q_i)} w(q_i, q_j)$$

含义：该 qubit 在 DD 中向上/向下各有多少"引力"。引力大的方向应优先 sift。

---

## 4. 变量选择策略

### 策略 A：纯中心性排序

按加权度降序处理变量。交互最多的 qubit 对全局 DD size 影响最大，先调整它们。

### 策略 B：混合评分（节点数 × 中心性）

$$\text{Score}(q_i) = \alpha \cdot \frac{\text{active}[\text{varMap}[i]]}{\max_j \text{active}[j]} + (1-\alpha) \cdot \frac{C_D(q_i)}{\max_j C_D(q_j)}$$

用 $\alpha = 0.6$ 兼顾当前 DD 状态和电路结构先验。

### 策略 C：引力决定 sifting 方向

替代原来 `pos < n/2 → 先下后上` 的启发式：

- $\text{Gravity}_\text{up} > \text{Gravity}_\text{down}$ → 先向上 sift
- 否则 → 先向下 sift

让变量优先向交互伙伴密集的方向移动。

---

## 5. 理论依据

让高权边对应的 qubit 在 DD 中尽量相邻，等价于经典的 **Minimum Linear Arrangement (MLA)** 问题：

$$\min_\pi \sum_{(u,v) \in E} w(u,v) \cdot |\pi(u) - \pi(v)|$$

MLA 是 NP-hard，但 interaction graph 驱动的 sifting 通过贪心策略逼近最优解。

---

## 6. 复杂度分析

| 操作 | 复杂度 |
|------|--------|
| 构建 interaction graph | $O(|gates| \cdot k^2)$，$k$ = 每门 qubit 数 |
| 计算加权度 | $O(n^2)$ |
| sifting 顺序排序 | $O(n \log n)$ |
| 方向引力计算 | $O(n)$ per variable |
| 总额外开销 | 相比 sifting 本体 $O(n^2 \cdot \text{exchange\_cost})$ 可忽略 |

---

## 7. 与现有策略对比

| 策略 | 变量选择依据 | 方向决策 | 特点 |
|------|------------|----------|------|
| 标准 Sifting | max active 节点数 | pos < n/2 先下后上 | 无先验 |
| IG-Degree | 加权度中心性 | 引力方向 | 利用电路结构 |
| IG-混合 | 节点数 × 中心性 | 引力方向 | 兼顾两者 |

---

## 8. 代码集成

### 新增文件

- `extern/dd_package/include/InteractionGraph.h` — InteractionGraph 类定义
- `extern/dd_package/src/DDigSifting.cpp` — igSifting 实现

### 接口

```cpp
// DDpackage.h 中新增枚举值
enum DynamicReorderingStrategy {
    ...,
    IGSifting        // Interaction Graph 驱动的 Sifting
};

// DDpackage.h 中新增方法声明
std::tuple<Edge, unsigned int, unsigned int> igSifting(
    Edge in, std::map<unsigned short, unsigned short>& varMap,
    const InteractionGraph& ig);
```

### 使用示例

```cpp
qc::QuantumComputation qc;
qc.import("circuit.qasm");

auto dd = std::make_unique<dd::Package>();
auto varMap = qc.buildVarMap();

// 构建 interaction graph
dd::InteractionGraph ig;
ig.build(qc);

// 用 IG 驱动的 sifting
auto [result, minSize, maxSize] = dd->igSifting(func, varMap, ig);
```
