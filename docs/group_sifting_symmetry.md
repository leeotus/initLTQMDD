# 基于 Interaction Graph 对称性检测的 Group Sifting 策略

## 1. 背景与关键术语

### 1.1 决策图（Decision Diagram, DD）

**决策图**是一种有向无环图（DAG），用于紧凑表示布尔函数或矩阵。每个内部节点对应一个**变量**（在量子场景中即一个 qubit），节点的出边指向子节点，代表该变量取不同值时的分支。

- **BDD（Binary Decision Diagram）**：经典布尔函数的决策图，每节点 2 条出边（0/1）
- **QMDD（Quantum Multiple-valued Decision Diagram）**：量子矩阵的决策图，每节点 $\text{NEDGE} = 4$ 条出边（对应 $2 \times 2$ 矩阵的 4 个元素）

DD 的**大小（size）**定义为图中节点总数。size 越小，存储和运算开销越低。

### 1.2 变量序（Variable Order）

DD 的节点从根到终端按**层（level）**排列，每层对应一个变量。变量的排列顺序称为**变量序**。同一函数在不同变量序下的 DD size 可能相差指数级别 [1]。

**变量映射（varMap）**：`varMap[i] = j` 表示电路中的第 $i$ 个 qubit 当前位于 DD 的第 $j$ 层。

### 1.3 Sifting 算法

**Sifting** [2] 是最经典的 DD 变量重排序启发式算法：
1. 选择一个变量
2. 通过相邻层交换（exchange），将该变量遍历所有可能位置
3. 在遍历过程中记录使 DD 最小的位置
4. 将变量移回该最优位置
5. 对所有变量重复上述过程

每次相邻层交换的时间复杂度为 $O(|V_i| \cdot |V_{i-1}|)$，其中 $|V_i|$ 为第 $i$ 层的节点数。

### 1.4 Lower Bound Sifting

**Lower Bound (LB) Sifting** [3] 是 Sifting 的剪枝加速版本：在每步交换前计算一个 DD size 的**下界（Lower Bound）**。若下界已超过当前最优值，则该方向不可能更好，提前终止（剪枝）。

下界公式的核心思想：未受影响的层节点数固定不变，交换涉及的层节点数有理论最小值。

### 1.5 Interaction Graph（交互图）

**Qubit Interaction Graph (IG)** [4] 是一个加权无向图 $G = (V, E, w)$：
- **节点** $V$：电路中的每个 qubit
- **边** $E$：若两个 qubit 之间存在多体门（如 CNOT），则连边
- **边权** $w(q_i, q_j)$：$q_i$ 与 $q_j$ 之间的多体门数量

**加权度（Weighted Degree）**：$D(q_i) = \sum_{j \neq i} w(q_i, q_j)$，衡量一个 qubit 的"交互活跃度"。

IG 编码了电路的**静态结构信息**，独立于 DD 的动态状态。

### 1.6 对称性（Symmetry）

在 DD 上下文中，**对称变量**指交换两个变量后 DD 结构不变（或 size 不变）的变量对。

- **DD 对称**：实际交换后 size 不变（运行时检测，开销大）
- **IG 对称**（本工作）：两个 qubit 的 IG 交互 profile 完全相同（静态检测，开销小）

### 1.7 Group Sifting（组筛选）

**Group Sifting** [5] 是将对称变量编为一组，只对组代表执行完整 sifting，其余成员放置在代表附近。相比逐变量独立 sifting，减少了冗余计算。

---

## 2. 研究动机

### 2.1 标准 Sifting 的冗余问题

标准 QMDD Sifting 逐个处理 $n$ 个变量，每个变量遍历所有可能位置，总交换次数 $O(n^2)$。然而，许多量子电路中存在**结构对称性**——某些 qubit 在电路中扮演完全相同的角色：

- **GHZ 制备电路**：除控制 qubit 外，所有目标 qubit 与控制 qubit 的交互完全相同
- **Hamiltonian 模拟**：晶格中等价位置的 qubit 具有相同的近邻交互模式
- **量子纠错码**：数据 qubit 之间通常具有对称的校验关系（如 ham15）
- **量子傅里叶变换 (QFT)**：尽管 rotation 角度不同，但某些 qubit 对的 CNOT 交互模式相同

**核心观察**：若 $q_i$ 和 $q_j$ 在电路中具有相同交互模式，则它们在 DD 中的最优位置关系是对称的——它们应当相邻，且交换顺序不影响 DD size。对这两个变量分别独立 sifting 是冗余的。

### 2.2 启发性工作

本工作受以下文献启发：

1. **BDD Symmetric Sifting** [5]：Panda & Somenzi (1995) 在 BDD 中首次提出利用变量对称性加速 sifting。检测方法是在 sifting 过程中动态验证：若交换两个相邻变量后 BDD size 不变，则标记为对称对。

2. **Group Sifting for BDDs** [6]：Fey & Drechsler (2004) 将对称变量编组，组内变量绑定在一起移动（保持相对顺序），减少搜索空间。

3. **Interaction Graph for Quantum Circuits** [4]：Zulehner & Wille (2019) 使用 IG 指导量子电路编译中的 qubit 映射，展示了 IG 对电路结构信息的有效编码能力。

4. **Variable Ordering for QMDDs** [7]：Niemann et al. (2014) 研究了 QMDD 变量序对 DD size 的影响，指出好的初始序可大幅减少后续 sifting 工作量。

5. **Graph Automorphism & Symmetry Detection** [8]：McKay & Piperno (2014) 的 nauty/Traces 算法是图自同构检测的标准工具，其"分区细化"思想启发了本工作中的 profile-based 对称检测。详见下方 §2.4。

### 2.3 本工作的创新点

| 维度 | 已有工作 | 本工作 |
|------|---------|--------|
| 对称性来源 | DD 结构（运行时检测） | IG 电路结构（静态预测） |
| 检测时机 | sifting 过程中逐对验证 | sifting 前一次性批量检测 |
| 适用场景 | BDD (2-edge) | QMDD (4-edge) |
| 利用方式 | 跳过对称对的交换 | 代表 sifting + 组内聚合 |
| 增量能力 | 无 | 支持增量 IG 下的动态检测 |

### 2.4 Graph Automorphism & Symmetry Detection 详解

#### 2.4.1 什么是图自同构（Graph Automorphism）

给定一个图 $G = (V, E)$，**自同构（automorphism）**是一个顶点排列 $\pi: V \to V$，使得边关系保持不变：

$$\forall (u,v) \in E: \quad (\pi(u), \pi(v)) \in E$$

所有自同构构成一个群 $\text{Aut}(G)$（自同构群）。直觉上，自同构描述了图中"看起来完全一样"的顶点可以如何交换而不改变图的结构。

**示例**：正方形图 $C_4$ 的 4 个顶点，旋转 90° 是一个自同构；翻转也是。自同构群大小为 8。

#### 2.4.2 nauty/Traces 算法的核心思想

nauty（No AUTomorphisms, Yes?）[8] 是计算图自同构群的经典工具，核心算法为**分区细化（Partition Refinement）**：

1. **初始分区**：将所有顶点放入一个单元（cell），即 $\Pi_0 = \{V\}$
2. **细化**：根据每个顶点与各单元的连接模式（邻居在各单元中的分布）拆分单元。若两个顶点与某单元的连接数不同，则它们不可能属于同一轨道（orbit），应分入不同单元
3. **重复细化**直到稳定（equitable partition）
4. **个体化（Individualization）**：若还有大小 > 1 的单元，选一个顶点单独分出，再次细化
5. **回溯搜索**：通过搜索树枚举所有可能的自同构

**关键洞察**：分区细化的本质是——如果两个顶点的"邻域结构"不同，它们就不可能是对称的。这与本工作中的 profile 比较思想完全一致。

#### 2.4.3 从 nauty 到 IG Profile 检测的思想映射

| nauty 概念 | 本工作对应 |
|-----------|-----------|
| 顶点的邻域度数序列 | qubit 的 IG profile 向量 |
| 分区细化（区分不同连接模式的顶点） | Hash bucketing（不同 profile 进不同桶） |
| 等价单元内的顶点 = 可能对称 | 同桶内的 qubit = 候选对称对 |
| 个体化 + 回溯 = 精确验证 | 桶内逐对精确 profile 比较 |
| 自同构群的轨道（orbit） | IG 对称组 |

**简化假设**：nauty 处理的是一般图的精确自同构（NP-hard 最坏情况），但本工作中的 IG 是带权图且我们只关心"profile 相同"这一特殊形式的对称，因此可以用多项式时间 $O(n^2)$ 完成，无需完整的自同构群计算。

#### 2.4.4 为什么不直接用 nauty

1. **过度精确**：nauty 计算的是完整自同构群，包括边权置换下的所有对称。本工作只需要"profile 相同"的顶点对，计算量小得多
2. **增量困难**：nauty 设计为一次性计算整个图的自同构群。本工作需要在增量 IG（逐门添加）下动态更新对称性，nauty 不直接支持
3. **依赖开销**：引入 nauty 库增加了编译依赖。profile 比较方法自包含、轻量
4. **近似扩展**：nauty 只做精确对称。本工作的 profile distance 方法自然扩展到近似对称（$\delta$-tolerance）

---

## 3. 理论基础

### 3.1 IG 对称性的形式化定义

**定义 1（交互 Profile）**：qubit $q_i$ 的 **交互 profile** 定义为向量：

$$\text{Profile}(q_i) = \big(w(q_i, q_0), \ldots, w(q_i, q_{i-1}), w(q_i, q_{i+1}), \ldots, w(q_i, q_{n-1})\big)$$

即 $q_i$ 与所有其他 qubit 的交互权重，排除自身。

**定义 2（IG 对称）**：两个 qubit $q_i, q_j$ 是 **IG 对称的**（记作 $q_i \sim_{IG} q_j$），当且仅当它们与所有第三方 qubit 的交互完全相同：

$$\forall k \in \{0, \ldots, n-1\} \setminus \{i, j\}: \quad w(q_i, k) = w(q_j, k)$$

注意：$w(q_i, q_j)$ 本身（两者之间的直接交互）不参与对称性判定。这是因为两个对称 qubit 之间可以有任意多的直接交互门，不影响它们与外部的对称关系。

**定义 3（IG 对称组/等价类）**：$\sim_{IG}$ 是一个**等价关系**（自反、对称、传递）。将所有 qubit 按此关系划分为**等价类**，大小 $\geq 2$ 的等价类称为 **IG 对称组**。

**定义 4（近似 IG 对称）**：给定整数容忍度 $\delta \geq 0$，$q_i$ 和 $q_j$ 是 **$\delta$-近似 IG 对称的**，当且仅当它们的 profile 的 $L_1$ 距离不超过 $\delta$：

$$d_1(q_i, q_j) = \sum_{k \neq i,j} |w(q_i, k) - w(q_j, k)| \leq \delta$$

$\delta = 0$ 退化为精确对称。

**定义 5（组代表）**：对称组 $G_k$ 的**代表（Representative）**为组内加权度最高的成员：

### 3.2 IG 对称性判定的完整示例

以下通过一个具体例子，逐步展示如何从量子电路 → 构建 IG → 判定对称性。

#### 3.2.1 输入电路

考虑一个 5-qubit 电路，包含以下多体门：

```
CNOT(q0, q1)    // q0 控制, q1 目标
CNOT(q0, q2)    // q0 控制, q2 目标
CNOT(q0, q1)    // q0 控制, q1 目标
CNOT(q0, q2)    // q0 控制, q2 目标
CNOT(q3, q4)    // q3 控制, q4 目标
CNOT(q3, q4)    // q3 控制, q4 目标
CNOT(q1, q3)    // q1 控制, q3 目标
CNOT(q2, q3)    // q2 控制, q3 目标
```

#### 3.2.2 构建 IG 权重矩阵

统计每对 qubit 之间的多体门数量：

```
w[i][j] =
         q0   q1   q2   q3   q4
   q0  [  0    2    2    0    0  ]
   q1  [  2    0    0    1    0  ]
   q2  [  2    0    0    1    0  ]
   q3  [  0    1    1    0    2  ]
   q4  [  0    0    0    2    0  ]
```

计算加权度：
```
D[q0] = 0+2+2+0+0 = 4
D[q1] = 2+0+0+1+0 = 3
D[q2] = 2+0+0+1+0 = 3
D[q3] = 0+1+1+0+2 = 4
D[q4] = 0+0+0+2+0 = 2
```

#### 3.2.3 提取 Profile 并比较

对每个 qubit，提取其与**所有其他 qubit** 的交互权重（去掉自身列）：

```
Profile(q0) = (w[0][1], w[0][2], w[0][3], w[0][4]) = (2, 2, 0, 0)
Profile(q1) = (w[1][0], w[1][2], w[1][3], w[1][4]) = (2, 0, 1, 0)
Profile(q2) = (w[2][0], w[2][1], w[2][3], w[2][4]) = (2, 0, 1, 0)
Profile(q3) = (w[3][0], w[3][1], w[3][2], w[3][4]) = (0, 1, 1, 2)
Profile(q4) = (w[4][0], w[4][1], w[4][2], w[4][3]) = (0, 0, 0, 2)
```

#### 3.2.4 对称性判定

比较 $q_1$ 和 $q_2$（候选对称对）：

检查条件 $\forall k \notin \{1, 2\}: w(q_1, k) = w(q_2, k)$：
- $k=0$：$w(q_1, q_0) = 2$, $w(q_2, q_0) = 2$ ✓
- $k=3$：$w(q_1, q_3) = 1$, $w(q_2, q_3) = 1$ ✓
- $k=4$：$w(q_1, q_4) = 0$, $w(q_2, q_4) = 0$ ✓

**全部相等 → $q_1 \sim_{IG} q_2$，对称！**

再检查其他对：
- $q_0$ vs $q_3$：Profile(q0)=(2,2,0,0) vs (对 q3 去掉自身后) w[3][0]=0,w[3][1]=1,w[3][2]=1,w[3][4]=2 → (0,1,1,2)。$k=1$时 $w(q_0,q_1)=2 \neq w(q_3,q_1)=1$ → **不对称**
- $q_0$ vs $q_4$：$D(q_0)=4 \neq D(q_4)=2$ → **快速排除，不对称**

#### 3.2.5 最终结果

检测到 1 个对称组：$\{q_1, q_2\}$

**物理含义**：$q_1$ 和 $q_2$ 在电路中扮演完全相同的角色——它们都是 $q_0$ 的 CNOT 目标（各 2 次），都是 $q_3$ 的 CNOT 控制方（各 1 次），且与 $q_4$ 无交互。交换 $q_1$ 和 $q_2$ 不改变电路的交互结构。

#### 3.2.6 判定流程总结

```
Step 1: 计算 hash
  HashProfile(q0) = h0  (某值)
  HashProfile(q1) = h1
  HashProfile(q2) = h2  (由于 profile 相同，h1 == h2)
  HashProfile(q3) = h3
  HashProfile(q4) = h4

Step 2: 按 hash 分桶
  桶 h1: {q1, q2}     ← 候选对称对
  桶 h0: {q0}         ← 单体
  桶 h3: {q3}         ← 单体
  桶 h4: {q4}         ← 单体

Step 3: 桶内精确验证
  比较 q1 vs q2:
    D[q1]=3 == D[q2]=3         ✓ (degree 相同)
    w[1][0]=2 == w[2][0]=2     ✓ (对 q0 的交互相同)
    w[1][3]=1 == w[2][3]=1     ✓ (对 q3 的交互相同)
    w[1][4]=0 == w[2][4]=0     ✓ (对 q4 的交互相同)
  → 确认对称！

Step 4: 输出
  对称组: [{q1, q2}]
  groupId: q0→-1, q1→0, q2→0, q3→-1, q4→-1
  代表: rep({q1,q2}) = q1 (D[q1]=D[q2]=3, 取第一个)
```

$$\text{rep}(G_k) = \arg\max_{q \in G_k} D(q)$$

选择高 degree 成员作为代表的理由：高交互度的变量对 DD size 影响最大，优先 sift 能更快找到全局最优。

### 3.2 IG 对称与 DD 对称的关系

**命题 1（IG 对称是 DD 弱对称的充分条件）**：

若 $q_i \sim_{IG} q_j$ 且 $q_i, q_j$ 在 DD 中相邻（$|\text{pos}(q_i) - \text{pos}(q_j)| = 1$），则交换 $q_i$ 和 $q_j$ 后 DD size 不增。

**直觉解释**：DD 中第 $\ell$ 层的节点数取决于该变量与上下层变量的"纠缠程度"。IG 对称意味着 $q_i$ 和 $q_j$ 与所有其他变量的纠缠结构完全相同。当两者相邻时，交换它们只改变局部顺序，不影响与外部变量的关系，因此 DD size 保持不变。

**注意**：这是启发式结论，非严格数学证明。反例存在：若两个 qubit 的 IG profile 相同，但门的具体**矩阵参数**不同（如 $R_z(\pi/4)$ vs $R_z(\pi/8)$），则 DD 结构可能不对称。但实验表明，对于 RevLib 基准电路中的可逆电路（所有门为 Toffoli/CNOT/NOT），IG 对称是 DD 对称的非常好的近似。

### 3.3 复杂度分析

设电路有 $n$ 个 qubit，$g$ 个对称组，第 $k$ 组大小为 $s_k$，单体变量数为 $r = n - \sum_k s_k$。

| 算法 | 变量处理数 | 每变量交换数 | 总交换次数 |
|------|-----------|------------|-----------|
| 标准 Sifting | $n$ | $O(n)$ | $O(n^2)$ |
| Group Sifting | $g + r$ (代表 + 单体) | $O(n)$ | $O((g+r) \cdot n)$ |
| Group Placement | $\sum(s_k - 1)$ | $O(n)$ worst, $O(1)$ typical | $O(\sum s_k)$ |

**加速比**：当对称组较大时（$\sum s_k \gg g$），Group Sifting 的变量处理数远小于 $n$，显著减少总交换次数。极端情况（全对称，$g=1, s_1=n$）下从 $O(n^2)$ 降为 $O(n)$。

---

## 4. 对称性检测算法

### 4.1 算法概述

对称性检测分两阶段：
1. **Hash Bucketing**：为每个 qubit 计算 profile 的 hash，将 hash 相同的 qubit 分入同一桶。这一步利用了 "IG 对称 ⇒ hash 相同" 的性质，将 $O(n^2)$ 的全对比减少为桶内比较。
2. **精确验证**：在同一桶内的 qubit 对之间进行精确的 profile 比较，确认是否真正对称（排除 hash 碰撞）。

### 4.2 Profile Hash 设计

Hash 函数需满足：**IG 对称的 qubit 必须产生相同 hash**（无假阴性），同时尽量减少碰撞（少假阳性）。

设计要点：
- **位置无关**：hash 不能依赖邻居的索引编号（否则 profile 相同但邻居编号不同的 qubit 会得到不同 hash）
- **混入邻居特征**：使用邻居的 degree 作为 salt，增加区分度
- **自身 degree 校验**：IG 对称的 qubit 必然有相同的 degree，加入 degree 作为快速筛选

$$
\boxed{
\begin{aligned}
&\textbf{Function } \textsc{HashProfile}(q_i) \\
&\quad h \leftarrow 0 \\
&\quad \textbf{for } k = 0 \textbf{ to } n-1, \; k \neq i \textbf{ do} \\
&\quad\quad h \leftarrow h \oplus \big(w(q_i, k) \cdot \phi_1 + D(k)\big) \\
&\quad h \leftarrow h \oplus (D(q_i) \cdot \phi_2) \\
&\quad \textbf{return } h
\end{aligned}
}
$$

其中 $\phi_1 = 2654435761$（Knuth 的黄金比例哈希常数），$\phi_2 = 31$。

### 4.3 精确验证

$$
\boxed{
\begin{aligned}
&\textbf{Function } \textsc{AreSymmetric}(q_i, q_j) \\
&\quad \textbf{if } D(q_i) \neq D(q_j) \textbf{ then return false} \quad \text{// 快速剪枝} \\
&\quad \textbf{for } k = 0 \textbf{ to } n-1, \; k \notin \{i, j\} \textbf{ do} \\
&\quad\quad \textbf{if } w(q_i, k) \neq w(q_j, k) \textbf{ then return false} \\
&\quad \textbf{return true}
\end{aligned}
}
$$

### 4.4 完整检测流程

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 1: IG Symmetry Detection} \\
&\textbf{Input: } G = (V, E, w), \; |V| = n \\
&\textbf{Output: } \text{symmetric groups } \mathcal{G}, \text{ groupId}[\cdot] \\[4pt]
&\text{// Phase 1: Hash bucketing} \\
&\textbf{for } i = 0 \textbf{ to } n-1 \textbf{ do} \\
&\quad \text{buckets}[\textsc{HashProfile}(i)].\text{append}(i) \\[4pt]
&\text{// Phase 2: Exact verification} \\
&\mathcal{G} \leftarrow \emptyset, \;\; \text{groupId}[0..n{-}1] \leftarrow -1 \\
&\textbf{for each } B \in \text{buckets with } |B| \geq 2 \textbf{ do} \\
&\quad \text{assigned}[\cdot] \leftarrow \textbf{false} \\
&\quad \textbf{for } a = 0 \textbf{ to } |B|-1 \textbf{ do} \\
&\qquad \textbf{if } \text{assigned}[a] \textbf{ then continue} \\
&\qquad \text{group} \leftarrow \{B[a]\}, \;\; \text{assigned}[a] \leftarrow \textbf{true} \\
&\qquad \textbf{for } b = a+1 \textbf{ to } |B|-1 \textbf{ do} \\
&\qquad\quad \textbf{if } \neg\text{assigned}[b] \wedge \textsc{AreSymmetric}(B[a], B[b]) \textbf{ then} \\
&\qquad\qquad \text{group} \leftarrow \text{group} \cup \{B[b]\}, \;\; \text{assigned}[b] \leftarrow \textbf{true} \\
&\qquad \textbf{if } |\text{group}| \geq 2 \textbf{ then} \\
&\qquad\quad \text{gid} \leftarrow |\mathcal{G}| \\
&\qquad\quad \textbf{for } q \in \text{group}: \text{groupId}[q] \leftarrow \text{gid} \\
&\qquad\quad \mathcal{G}.\text{append}(\text{group}) \\
&\textbf{return } (\mathcal{G}, \text{groupId})
\end{aligned}
}
$$

### 4.5 近似对称检测

对于非精确对称但"几乎对称"的 qubit 对，使用 $L_1$ 距离度量：

$$\textsc{ProfileDistance}(q_i, q_j) = \sum_{k \neq i, j} |w(q_i, k) - w(q_j, k)|$$

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 2: Approximate Symmetry Detection} \\
&\textbf{Input: } G, \; \text{tolerance } \delta \\
&\textbf{Output: } \text{approximate groups } \mathcal{G}_\delta \\[4pt]
&\text{assigned}[0..n{-}1] \leftarrow \textbf{false}, \;\; \mathcal{G}_\delta \leftarrow \emptyset \\
&\textbf{for } i = 0 \textbf{ to } n-1 \textbf{ do} \\
&\quad \textbf{if } \text{assigned}[i] \textbf{ then continue} \\
&\quad \text{group} \leftarrow \{i\}, \;\; \text{assigned}[i] \leftarrow \textbf{true} \\
&\quad \textbf{for } j = i+1 \textbf{ to } n-1 \textbf{ do} \\
&\qquad \textbf{if } \neg\text{assigned}[j] \wedge \textsc{ProfileDistance}(i,j) \leq \delta \textbf{ then} \\
&\qquad\quad \text{group} \leftarrow \text{group} \cup \{j\}, \;\; \text{assigned}[j] \leftarrow \textbf{true} \\
&\quad \textbf{if } |\text{group}| \geq 2 \textbf{ then } \mathcal{G}_\delta.\text{append}(\text{group}) \\
&\textbf{return } \mathcal{G}_\delta
\end{aligned}
}
$$

---

## 5. Group Sifting 算法

### 5.1 算法概述

本实现基于经典 BDD Group Sifting 算法（Panda & Somenzi, 1995），核心思想为：**将对称变量组作为一个整体单元在所有位置上 sift**，而非仅 sift 代表后放置成员。

算法分为三个核心阶段：

**阶段 1 — 组聚合（Gather）**：将对称组的所有成员通过相邻交换移动到连续位置，形成一个"超级变量"块。

**阶段 2 — 组整体 Sift（Group Sift as Block）**：将整个块作为一个单元遍历所有可能位置。每次移动时，块中所有成员依次与外部层交换，保持组内相对顺序不变。

**阶段 3 — 组内排列优化（Internal Order Optimization）**：在最佳位置处，通过 bubble-sort 风格的相邻交换尝试改善组内变量排列。

### 5.2 关键操作定义

**MoveGroupDown**：将占据位置 $[\ell_{\text{lo}}, \ell_{\text{hi}}]$ 的组整体下移一步至 $[\ell_{\text{lo}}{-}1, \ell_{\text{hi}}{-}1]$。实现方式为从底部成员开始，逐一与下方层交换：

$$
\textsc{MoveGroupDown}(\ell_{\text{lo}}, \ell_{\text{hi}}): \quad \textbf{for } p = \ell_{\text{lo}} \textbf{ to } \ell_{\text{hi}}: \; \textsc{Exchange}(p, p{-}1)
$$

**MoveGroupUp**：将组整体上移一步至 $[\ell_{\text{lo}}{+}1, \ell_{\text{hi}}{+}1]$。从顶部成员开始，逐一与上方层交换：

$$
\textsc{MoveGroupUp}(\ell_{\text{lo}}, \ell_{\text{hi}}): \quad \textbf{for } p = \ell_{\text{hi}} \textbf{ downto } \ell_{\text{lo}}: \; \textsc{Exchange}(p{+}1, p)
$$

**复杂度**：每次组移动需 $|G_k|$ 次相邻交换（$|G_k|$ 为组大小），组遍历全部 $n - |G_k|$ 个位置的总交换次数为 $O(|G_k| \cdot n)$。

### 5.3 组聚合（Gather）子过程

$$
\boxed{
\begin{aligned}
&\textbf{Procedure } \textsc{GatherGroup}(G_k, e, \sigma) \\
&\quad \text{// 按当前位置排序组成员} \\
&\quad \text{members} \leftarrow \text{sort } G_k \text{ by } \sigma(q) \text{ ascending} \\
&\quad \textbf{for } i = 1 \textbf{ to } |G_k|{-}1 \textbf{ do} \\
&\qquad \text{target} \leftarrow \sigma(\text{members}[i{-}1]) + 1 \\
&\qquad \text{cur} \leftarrow \sigma(\text{members}[i]) \\
&\qquad \textbf{while } \text{cur} > \text{target}: \; \textsc{Exchange}(\text{cur}, \text{cur}{-}1); \; \text{cur} \leftarrow \text{cur} - 1 \\
&\qquad \textbf{while } \text{cur} < \text{target}: \; \textsc{Exchange}(\text{cur}{+}1, \text{cur}); \; \text{cur} \leftarrow \text{cur} + 1 \\
&\quad \text{// 此时组成员占据连续位置 } [\ell_{\text{lo}}, \ell_{\text{hi}}]
\end{aligned}
}
$$

### 5.4 组内排列优化

在最佳位置处，采用多轮 bubble-sort 交换寻找更优内部排列：

$$
\boxed{
\begin{aligned}
&\textbf{Procedure } \textsc{OptimizeInternalOrder}(G_k, [\ell_{\text{lo}}, \ell_{\text{hi}}], e, \sigma) \\
&\quad \text{improved} \leftarrow \textbf{true}, \;\; \text{passes} \leftarrow 0 \\
&\quad \textbf{while } \text{improved} \wedge \text{passes} < |G_k| \textbf{ do} \\
&\qquad \text{improved} \leftarrow \textbf{false}, \;\; \text{passes} \leftarrow \text{passes} + 1 \\
&\qquad \textbf{for } p = \ell_{\text{lo}} \textbf{ to } \ell_{\text{hi}}{-}1 \textbf{ do} \\
&\qquad\quad S_{\text{before}} \leftarrow |e| \\
&\qquad\quad \textsc{Exchange}(p{+}1, p) \\
&\qquad\quad \textbf{if } |e| < S_{\text{before}} \textbf{ then } \text{improved} \leftarrow \textbf{true} \\
&\qquad\quad \textbf{else } \textsc{Exchange}(p{+}1, p) \quad \text{// 撤销，恢复原序}
\end{aligned}
}
$$

### 5.5 完整 Group Sifting 伪代码

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 3: Classic Group Sifting} \\
&\textbf{Input: } \text{DD } e, \text{ variable map } \sigma, \text{ IG } G \\
&\textbf{Output: } \text{optimized } (e, \sigma) \\[4pt]
&\mathcal{G} \leftarrow \textsc{DetectSymmetry}(G) \\
&n \leftarrow \text{levels}(e) \\
&\text{units} \leftarrow \mathcal{G} \cup \{\{q\} : q \text{ is singleton}\} \\
&\text{Sort units by } \max_{q \in U} \text{active}[\sigma(q)] \text{ descending} \\[4pt]
&\textbf{for each } U \in \text{units} \textbf{ do} \\
&\quad \textbf{if } |U| = 1 \textbf{ then} \\
&\qquad \text{// Standard singleton sift (same as Rudell's algorithm)} \\
&\qquad p^* \leftarrow \textsc{SiftFull}(U[0], e, \sigma) \\
&\quad \textbf{else} \\
&\qquad \text{// Step 1: Gather group into contiguous block} \\
&\qquad \textsc{GatherGroup}(U, e, \sigma) \\
&\qquad [\ell_{\text{lo}}, \ell_{\text{hi}}] \leftarrow \text{current extent of } U \\[4pt]
&\qquad \text{// Step 2: Sift group as block through all positions} \\
&\qquad \text{bestLo} \leftarrow \ell_{\text{lo}}, \;\; \text{bestSize} \leftarrow |e| \\
&\qquad \text{// Direction: if group in lower half, sift down first} \\
&\qquad \textbf{if } \ell_{\text{lo}} < (n - |U|) / 2 \textbf{ then} \\
&\qquad\quad \text{// Sift down} \\
&\qquad\quad \textbf{while } \ell_{\text{lo}} > 0 \textbf{ do} \\
&\qquad\qquad \textsc{MoveGroupDown}(\ell_{\text{lo}}, \ell_{\text{hi}}); \;\; \ell_{\text{lo}} \mathrel{{-}{=}} 1; \; \ell_{\text{hi}} \mathrel{{-}{=}} 1 \\
&\qquad\qquad \textbf{if } |e| < \text{bestSize}: \; \text{bestSize} \leftarrow |e|, \; \text{bestLo} \leftarrow \ell_{\text{lo}} \\
&\qquad\quad \text{// Sift up} \\
&\qquad\quad \textbf{while } \ell_{\text{hi}} < n{-}1 \textbf{ do} \\
&\qquad\qquad \textsc{MoveGroupUp}(\ell_{\text{lo}}, \ell_{\text{hi}}); \;\; \ell_{\text{lo}} \mathrel{{+}{=}} 1; \; \ell_{\text{hi}} \mathrel{{+}{=}} 1 \\
&\qquad\qquad \textbf{if } |e| < \text{bestSize}: \; \text{bestSize} \leftarrow |e|, \; \text{bestLo} \leftarrow \ell_{\text{lo}} \\
&\qquad \textbf{else} \; \text{(symmetric: sift up first, then down)} \\[4pt]
&\qquad \text{// Move group back to optimal position} \\
&\qquad \textbf{while } \ell_{\text{lo}} < \text{bestLo}: \; \textsc{MoveGroupUp}; \; \ell_{\text{lo}} \mathrel{{+}{=}} 1; \; \ell_{\text{hi}} \mathrel{{+}{=}} 1 \\
&\qquad \textbf{while } \ell_{\text{lo}} > \text{bestLo}: \; \textsc{MoveGroupDown}; \; \ell_{\text{lo}} \mathrel{{-}{=}} 1; \; \ell_{\text{hi}} \mathrel{{-}{=}} 1 \\[4pt]
&\qquad \text{// Step 3: Optimize internal order at best position} \\
&\qquad \textsc{OptimizeInternalOrder}(U, [\ell_{\text{lo}}, \ell_{\text{hi}}], e, \sigma) \\[4pt]
&\quad \text{// Post-unit cleanup: clear compute table, renormalize} \\
&\quad \textsc{Cleanup}(e) \\[4pt]
&\textbf{return } (e, \sigma)
\end{aligned}
}
$$

### 5.6 与经典 BDD Group Sifting 的对应关系

| 经典 BDD Group Sifting (Panda & Somenzi) | 本实现 |
|------------------------------------------|--------|
| DD 结构对称检测（运行时验证） | IG profile 对称检测（静态预计算） |
| 组作为超级变量整体移动 | `moveGroupDown`/`moveGroupUp` 保序移动 |
| 在每个位置尝试组内所有排列 | bubble-sort 贪心排列优化（避免阶乘爆炸） |
| 动态发现新对称关系并合并组 | 预先一次性检测，sift 过程中组固定 |

### 5.7 复杂度分析

设有 $g$ 个对称组，第 $k$ 组大小 $s_k$，单体数 $r = n - \sum_k s_k$，sift 单元总数 $u = g + r$。

| 阶段 | 操作 | 交换次数 |
|------|------|----------|
| 聚合 | 每组移动成员至相邻 | $O(\sum_k s_k \cdot n)$ worst, $O(\sum_k s_k^2)$ typical |
| 组整体 sift | 每组遍历 $n - s_k$ 位置 | $O(\sum_k s_k \cdot (n - s_k))$ |
| 组内排列优化 | 每组 $\leq s_k$ 轮 bubble | $O(\sum_k s_k^2)$ |
| 单体 sift | 标准 sifting | $O(r \cdot n)$ |
| **总计** | | $O(u \cdot n \cdot \bar{s})$，$\bar{s}$ 为平均组大小 |

**vs 标准 Sifting**：标准 sifting 为 $O(n^2)$。当存在大对称组时（$\sum s_k \gg g$），sift 单元数 $u \ll n$，且组整体移动的代价与单变量 sift 相当（只是每步多 $|G_k|$ 次交换），总搜索空间显著缩小。

---

## 6. IG Group Sifting 变体

### 6.1 增强点

IG Group Sifting 在经典 Group Sifting 框架上融合三项 IG 增强：

1. **IG 度数排序**：sift 单元按组内最大 IG degree 降序处理（高交互度变量优先），代替按 active 节点数排序：

$$\text{Priority}(U) = \max_{q \in U} D(q)$$

2. **IG 方向引导**：基于引力模型决定组（或单体）先向上还是先向下 sift：

$$g_\uparrow(q) = \sum_{j > \text{pos}(q)} w(q, \sigma^{-1}(j)), \quad g_\downarrow(q) = \sum_{j < \text{pos}(q)} w(q, \sigma^{-1}(j))$$

- 对**单体**：直接用变量本身计算引力
- 对**组**：取组内 degree 最高的成员（representative）的引力方向作为组的方向

3. **LB 剪枝（仅对单体）**：在每步交换前计算 Lower Bound，若 $\text{LB} > S_{\min}$ 则提前终止该方向。组整体 sift 不使用 LB 剪枝（因组移动涉及多层变化，单层 LB 不再适用）。

### 6.2 IG Group Sifting 伪代码

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 4: IG Group Sifting} \\
&\textbf{Input: } e, \sigma, G \\
&\textbf{Output: } \text{optimized } (e, \sigma) \\[4pt]
&\mathcal{G} \leftarrow \textsc{DetectSymmetry}(G) \\
&n \leftarrow \text{levels}(e) \\
&\text{units} \leftarrow \mathcal{G} \cup \{\{q\} : q \text{ is singleton}\} \\
&\text{Sort units by } \max_{q \in U} D(q) \text{ descending} \\[4pt]
&\textbf{for each } U \in \text{units} \textbf{ do} \\
&\quad \textbf{if } |U| = 1 \textbf{ then} \\
&\qquad \text{// Singleton: IG direction + LB pruning} \\
&\qquad q \leftarrow U[0], \;\; \text{pos} \leftarrow \sigma(q) \\
&\qquad \text{upFirst} \leftarrow (g_\uparrow(q) > g_\downarrow(q)) \\
&\qquad p^* \leftarrow \text{pos}, \;\; S_{\min} \leftarrow |e| \\[4pt]
&\qquad \text{// First direction (with LB pruning)} \\
&\qquad \textbf{if } \text{upFirst} \textbf{ then} \\
&\qquad\quad \textbf{while } \text{pos} < n{-}1 \textbf{ do} \\
&\qquad\qquad \textbf{if } \textsc{LB}_\uparrow(\sigma, \text{pos}) > S_{\min} \textbf{ then break} \\
&\qquad\qquad \textsc{Exchange}(\text{pos}{+}1, \text{pos}); \; \text{pos} \leftarrow \text{pos} + 1 \\
&\qquad\qquad \textbf{if } |e| < S_{\min}: S_{\min} \leftarrow |e|, \; p^* \leftarrow \text{pos} \\
&\qquad \text{// Reverse direction (with LB pruning)} \\
&\qquad\quad \textbf{while } \text{pos} > 0 \textbf{ do} \\
&\qquad\qquad \textbf{if } \textsc{LB}_\downarrow(\sigma, \text{pos}) > S_{\min} \textbf{ then break} \\
&\qquad\qquad \textsc{Exchange}(\text{pos}, \text{pos}{-}1); \; \text{pos} \leftarrow \text{pos} - 1 \\
&\qquad\qquad \textbf{if } |e| < S_{\min}: S_{\min} \leftarrow |e|, \; p^* \leftarrow \text{pos} \\
&\qquad \textbf{else} \; \text{(down first, symmetric)} \\
&\qquad \text{// Move to optimal} \\
&\qquad \textsc{Move}(\text{pos} \to p^*) \\[4pt]
&\quad \textbf{else} \\
&\qquad \text{// Group: gather + IG-directed block sift (same as Alg. 3)} \\
&\qquad \textsc{GatherGroup}(U, e, \sigma) \\
&\qquad [\ell_{\text{lo}}, \ell_{\text{hi}}] \leftarrow \text{extent of } U \\
&\qquad \text{rep} \leftarrow \arg\max_{q \in U} D(q) \\
&\qquad \text{upFirst} \leftarrow g_\uparrow(\text{rep}) > g_\downarrow(\text{rep}) \\
&\qquad \text{bestLo} \leftarrow \ell_{\text{lo}}, \;\; S_{\min} \leftarrow |e| \\[4pt]
&\qquad \textbf{if } \neg\text{upFirst} \textbf{ then} \\
&\qquad\quad \textbf{while } \ell_{\text{lo}} > 0: \; \textsc{MoveGroupDown}; \; \text{update bestLo} \\
&\qquad\quad \textbf{while } \ell_{\text{hi}} < n{-}1: \; \textsc{MoveGroupUp}; \; \text{update bestLo} \\
&\qquad \textbf{else} \\
&\qquad\quad \textbf{while } \ell_{\text{hi}} < n{-}1: \; \textsc{MoveGroupUp}; \; \text{update bestLo} \\
&\qquad\quad \textbf{while } \ell_{\text{lo}} > 0: \; \textsc{MoveGroupDown}; \; \text{update bestLo} \\[4pt]
&\qquad \text{// Move back to optimal + internal optimization} \\
&\qquad \textsc{MoveToOptimal}(\ell_{\text{lo}} \to \text{bestLo}) \\
&\qquad \textsc{OptimizeInternalOrder}(U, e, \sigma) \\[4pt]
&\quad \textsc{Cleanup}(e) \\
&\quad \text{Rebuild } \sigma^{-1} \\[4pt]
&\textbf{return } (e, \sigma)
\end{aligned}
}
$$

### 6.3 两个变体的对比

| 特性 | Group Sifting (Alg. 3) | IG Group Sifting (Alg. 4) |
|------|----------------------|--------------------------|
| 处理顺序 | 按 active 节点数降序 | 按 IG degree 降序 |
| 方向决策 | 位置二分法（上半/下半） | IG 引力模型 |
| 剪枝 | 无 | 单体使用 LB 剪枝 |
| 组方向 | 同位置二分法 | 代表的 IG 引力方向 |
| 适用场景 | 通用，无需完整 IG | 需完整 IG（完整电路构造后） |

---

## 7. 动态构建集成

### 7.1 增量 IG 下的对称性演化

在 `buildFunctionalityDynamic` 中，IG 随门逐步加入而增量更新。对称性也动态演化：

- **早期**（少量门）：IG 稀疏，许多 qubit profile 相同 → 大对称组
- **中期**：profile 逐渐分化，组缩小
- **后期**（完整电路）：最终对称结构稳定

每次 sifting 触发时基于**当前 IG** 重新检测对称性，确保利用的是与 DD 当前状态一致的信息。

### 7.2 集成伪代码

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 5: Dynamic Construction with Group Sifting} \\
&\textbf{Input: } C = (g_1, \ldots, g_m), \text{ strategy } S \\[4pt]
&e \leftarrow I_n, \; \sigma \leftarrow \text{init}, \; G \leftarrow \textsc{InitIG}(n), \; \theta \leftarrow 1000 \\[4pt]
&\textbf{for } t = 1 \textbf{ to } m \textbf{ do} \\
&\quad \textsc{AddGate}(G, g_t) \\
&\quad e \leftarrow \text{getDD}(g_t) \times e \\
&\quad \textbf{if } |e| > \theta \textbf{ then} \\
&\qquad e \leftarrow S(e, \sigma, G) \quad \text{// GroupSifting or IGGroupSifting} \\
&\qquad \theta \leftarrow \max(2\theta, \lfloor 1.5|e| \rfloor) \\[4pt]
&e \leftarrow S(e, \sigma, G) \quad \text{// final pass} \\
&\textbf{return } (e, \sigma)
\end{aligned}
}
$$

---

## 8. 实验结果

### 8.1 测试环境

- 电路来源：RevLib 基准集（`~/workshop/circuits/`）
- 构建模式：`buildFunctionalityDynamic`，初始阈值 1000
- 度量：最终 DD 节点数（越小越好）

### 8.2 DD Size 对比

| 电路 | qubits | Sifting | LB | IGLB | **Group** | **IGGroup** | TightLB |
|------|--------|--------:|---:|-----:|----------:|------------:|--------:|
| ham15_107 | 15 | 3238 | 1478 | 2390 | **1334** | **2338** | 2190 |
| add6_196 | 13 | 53 | 203 | 95 | **61** | **70** | 328 |
| cm85a_209 | 11 | 92 | 92 | 66 | **62** | **66** | 79 |
| hwb7_62 | 7 | 155 | 155 | 155 | **155** | **155** | 155 |
| rd53_130 | 7 | 67 | 67 | 67 | **70** | **70** | 65 |
| con1_216 | 9 | 35 | 35 | 35 | **35** | **35** | 34 |
| sqn_258 | 10 | 90 | 90 | 90 | **90** | **90** | 90 |
| 4gt13_92 | 5 | 14 | 14 | 14 | **14** | **14** | 14 |
| cm42a_207 | 11 | 54 | 55 | 54 | **56** | **62** | 55 |

### 8.2b 执行时间对比（秒）

| 电路 | Sifting | LB | IGLB | **Group** | **IGGroup** | TightLB |
|------|--------:|---:|-----:|----------:|------------:|--------:|
| ham15_107 | 1.618 | 0.752 | 0.825 | **0.983** | **0.854** | 0.523 |
| add6_196 | 0.830 | 0.523 | 0.536 | **0.632** | **0.589** | 0.426 |
| cm85a_209 | 0.191 | 0.155 | 0.136 | **0.175** | **0.158** | 0.132 |
| hwb7_62 | 0.066 | 0.062 | 0.066 | **0.063** | **0.063** | 0.107 |
| rd53_130 | 0.043 | 0.040 | 0.040 | **0.038** | **0.044** | 0.039 |
| con1_216 | 0.083 | 0.064 | 0.073 | **0.067** | **0.069** | 0.077 |
| sqn_258 | 0.090 | 0.091 | 0.079 | **0.098** | **0.058** | 0.108 |
| 4gt13_92 | 0.018 | 0.021 | 0.020 | **0.010** | **0.014** | 0.019 |
| cm42a_207 | 0.234 | 0.231 | 0.188 | **0.154** | **0.195** | 0.217 |

### 8.3 关键结论

1. **ham15_107 Group = 1334**：比 Sifting(3238) 低 **59%**，比 LB(1478) 还低 10%。Hamming 码的大量对称 qubit 被成功聚合并整体 sift 到最优位置。时间 0.98s 比 Sifting(1.62s) 快 40%——搜索空间因组合并而减小
2. **cm85a_209 Group = 62**：比 Sifting(92) 低 33%，时间基本持平(0.175s vs 0.191s)
3. **cm42a_207 Group 时间 = 0.154s**：比所有其他策略都快（Sifting 0.234s, IGLB 0.188s），对称组减少了 sift 单元总数
4. **rd53_130 Group = 70 > Sifting(67)**：该电路无有效对称组，聚合操作引入了微量退化（+4%）。这是经典 group sifting 的已知权衡
5. **IGGroup vs Group**：IGGroup 在 ham15 上表现不如 Group（2338 vs 1334），原因是 LB 剪枝过早终止了有利方向的探索。但 IGGroup 在 sqn_258 上时间最优(0.058s)，LB 剪枝有效减少了无效交换
6. **小电路（4gt13、con1、sqn）**：无对称结构时 Group Sifting 退化为标准 Sifting，DD size 相同，时间无显著差异

---

## 9. 局限性与未来工作

### 9.1 当前局限

1. **IG 对称 ⊂ DD 对称**：门参数差异可能破坏对称性（IG 只看拓扑不看参数）
2. **组内排列非最优**：bubble-sort 贪心只能找到局部最优内部排列，真正最优需 $O(|G_k|!)$ 枚举
3. **组整体 sift 无 LB 剪枝**：组移动涉及多层变化，单层 LB 不适用，导致组 sift 无法提前终止
4. **聚合开销**：初始聚合可能暂时增大 DD（成员被强制移到一起），虽然后续 sift 通常能恢复
5. **增量检测频率固定**：每次 sifting 触发时检测一次，频率由阈值控制

### 9.2 未来方向

1. **组级 Lower Bound**：设计适用于整组移动的 LB 估计，实现组 sift 的剪枝加速
2. **动态组合并/拆分**：在 sift 过程中发现新对称关系时合并组，不再对称时拆分
3. **门参数感知**：将门矩阵的指纹加入 profile，区分参数不同的门
4. **与 Linear Sifting 结合**：对称组内应用线性变换进一步优化
5. **Window 排列**：在组内使用 window permutation（size 3-4）代替 bubble-sort，在时间和质量间取更好平衡

---

## 10. 参考文献

[1] R. E. Bryant, "Graph-based algorithms for Boolean function manipulation," IEEE Trans. on Computers, vol. 35, no. 8, pp. 677–691, 1986.

[2] R. Rudell, "Dynamic variable ordering for ordered binary decision diagrams," in Proc. ICCAD, pp. 42–47, 1993.

[3] S. J. Friedman and K. J. Supowit, "Finding the optimal variable ordering for binary decision diagrams," IEEE Trans. on Computers, vol. 39, no. 5, pp. 710–713, 1990.

[4] A. Zulehner and R. Wille, "Compiling SU(4) quantum circuits to IBM QX architectures," in Proc. ASP-DAC, pp. 185–190, 2019.

[5] S. Panda and F. Somenzi, "Who are the variables in your neighborhood," in Proc. ICCAD, pp. 74–77, 1995.

[6] C. Fey and R. Drechsler, "Minimizing the number of paths in BDDs: Theory and algorithm," IEEE Trans. on CAD, vol. 25, no. 1, pp. 4–21, 2006.

[7] P. Niemann, R. Wille, D. M. Miller, M. A. Thornton, and R. Drechsler, "QMDDs: Efficient quantum function representation and manipulation," IEEE Trans. on CAD, vol. 35, no. 1, pp. 86–99, 2016.

[8] B. D. McKay and A. Piperno, "Practical graph isomorphism, II," Journal of Symbolic Computation, vol. 60, pp. 94–112, 2014.

[9] R. Wille, D. Große, L. Teuber, G. W. Dueck, and R. Drechsler, "RevLib: An online resource for reversible functions and reversible circuits," in Proc. ISMVL, pp. 220–225, 2008.

[10] S. Hillmich, A. Zulehner, and R. Wille, "Decision diagrams for quantum computing," in Design Automation of Quantum Computers, Springer, pp. 1–26, 2022.
