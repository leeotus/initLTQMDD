# Interaction Graph 通俗讲解

## 类比：社交关系网

把量子电路想象成一个**公司**，每个 qubit 是一个**员工**，每个双 qubit 门（如 CNOT）是两个员工之间的一次**合作**。

```
量子电路                        公司
─────────                      ────
qubit  q0, q1, q2, q3    →    员工 A, B, C, D
CNOT(q0, q1)              →    A 和 B 合作了 1 次
CNOT(q0, q1)              →    A 和 B 又合作了 1 次
CNOT(q1, q2)              →    B 和 C 合作了 1 次
CNOT(q0, q3)              →    A 和 D 合作了 1 次
```

**Interaction Graph 就是这个公司的"协作关系图"**——谁跟谁合作过，合作了几次。

---

## 用一个具体例子走一遍

假设有这样一个 4 qubit 电路：

```
q0: ──■──────■──────■──
      │      │      │
q1: ──X──■───X──────┼──
         │         │
q2: ─────X──────■──┼──
                │  │
q3: ────────────X──X──
```

门列表：
1. CNOT(q0, q1)
2. CNOT(q1, q2)
3. CNOT(q0, q1)
4. CNOT(q2, q3)
5. CNOT(q0, q3)

### 第一步：建图 —— 统计谁和谁交互了几次

```
           q0
          / | \
       2/   |  \1
       /    |   \
     q1     |   q3
       \    |   /
       1\   |  /1
         \  | /
          q2
```

写成矩阵：

```
weight[][] =
          q0   q1   q2   q3
    q0  [  0    2    0    1  ]     ← q0 和 q1 合作 2 次，和 q3 合作 1 次
    q1  [  2    0    1    0  ]     ← q1 和 q0 合作 2 次，和 q2 合作 1 次
    q2  [  0    1    0    1  ]     ← q2 和 q1 合作 1 次，和 q3 合作 1 次
    q3  [  1    0    1    0  ]     ← q3 和 q0 合作 1 次，和 q2 合作 1 次
```

### 第二步：算加权度 —— 每个 qubit 的"社交活跃度"

```
degree[q0] = 2 + 0 + 1 = 3    ← 最活跃！
degree[q1] = 2 + 1 + 0 = 3    ← 并列最活跃
degree[q2] = 0 + 1 + 1 = 2
degree[q3] = 1 + 0 + 1 = 2
```

### 第三步：用在 sifting 中

#### 变量选择

标准 sifting 选 active 节点最多的变量。IG 加入 degree 作为加分项：

```
Score(qi) = 0.85 × 节点数归一化 + 0.15 × degree归一化

假设 active 节点数: q0=10, q1=8, q2=12, q3=5

标准 sifting 选: q2 (active=12 最大)

IG sifting 计算:
  Score(q0) = 0.85×(10/12) + 0.15×(3/3) = 0.71 + 0.15 = 0.86
  Score(q1) = 0.85×(8/12)  + 0.15×(3/3) = 0.57 + 0.15 = 0.72
  Score(q2) = 0.85×(12/12) + 0.15×(2/3) = 0.85 + 0.10 = 0.95  ← 还是 q2
  Score(q3) = 0.85×(5/12)  + 0.15×(2/3) = 0.35 + 0.10 = 0.45

这里 q2 仍然胜出。但如果 q0 和 q2 的 active 接近时，IG 会让交互更多的 q0 优先。
```

#### 方向决策

假设 q0 当前在 DD 的 level 1（中间偏下）：

```
DD 层:   level 3 (顶)    level 2       level 1      level 0 (底)
变量:        q2              q3            q0             q1

q0 的交互伙伴在哪？
  q1 在 level 0 (下方)  → weight = 2
  q3 在 level 2 (上方)  → weight = 1

gravityDown = weight(q0, q1) = 2   ← q0 和 q1 合作频繁，q1 在下面
gravityUp   = weight(q0, q3) = 1

gravityDown > gravityUp → 先向下 sift！
```

**直觉**：把 q0 往下移、靠近 q1——因为它俩交互最多，**放相邻层可以最大化节点共享**。

---

## 为什么节点共享和"相邻"有关？

```
情况 A：q0 和 q1 相邻层              情况 B：q0 和 q1 隔了两层

┌────────────────────┐            ┌────────────────────┐
│ level 3: q2        │            │ level 3: q2        │
│ level 2: q3        │            │ level 2: q0        │  ← q0 和 q1 中间
│ level 1: q0        │ ← 相邻！   │ level 1: q3        │    隔了 q3
│ level 0: q1        │            │ level 0: q1        │
└────────────────────┘            └────────────────────┘

情况 A 中，q0 的子节点直接就是 q1 层的节点。
如果 CNOT(q0,q1) 使得 q0 的多个节点指向相同的 q1 子节点，
唯一表会自动合并 → DD 更小。

情况 B 中，q0 → q3 → q1，中间隔了 q3。
q0 的不同分支经过 q3 层时可能分裂出不同路径，
即使最终指向相同 q1 节点，中间的 q3 层节点也无法合并 → DD 更大。
```

**核心原理**：交互频繁的 qubit 放相邻层 = 更多共享 = DD 更小。

---

## 对应代码走查

### 构建 IG（`InteractionGraph.h`）

```cpp
template<typename QuantumCircuit>
void build(const QuantumCircuit& qc) {
    n = qc.getNqubits();                    // 有几个员工
    weight.assign(n, vector<int>(n, 0));    // 初始化合作次数矩阵

    for (const auto& op : qc) {             // 遍历每个门
        vector<unsigned short> involved;
        for (auto t : op->getTargets())     // 收集目标 qubit
            involved.push_back(t);
        for (auto& c : op->getControls())   // 收集控制 qubit
            involved.push_back(c.qubit);

        if (involved.size() < 2) continue;  // 单 qubit 门不产生交互

        // 对该门涉及的每一对 qubit，合作次数 +1
        for (size_t a = 0; a < involved.size(); ++a)
            for (size_t b = a + 1; b < involved.size(); ++b) {
                weight[involved[a]][involved[b]]++;
                weight[involved[b]][involved[a]]++;
            }
    }

    // 算每个 qubit 的总交互度
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            degree[i] += weight[i][j];
}
```

### 引力方向（`DDigSifting.cpp`）

```cpp
// 计算 pos 位置的变量向上/向下的"引力"
double gravityUp = 0, gravityDown = 0;

for (short j = pos + 1; j < n; ++j) {           // 遍历上方所有层
    auto circuitVar = invVarMap[j];               // 那一层放的是哪个 qubit？
    gravityUp += ig.weight[curVar][circuitVar];   // 和当前变量的交互次数
}

for (short j = 0; j < pos; ++j) {                // 遍历下方所有层
    auto circuitVar = invVarMap[j];
    gravityDown += ig.weight[curVar][circuitVar];
}

// 哪边引力大就先往哪边 sift
siftDownFirst = (gravityDown > gravityUp);
```

---

## 常见误区

**误区：IG 能让无剪枝 sifting 找到更好结果。**

错。无剪枝 sifting 每个变量遍历全部位置（从最底到最顶再回到最优点），不管先下还是先上，最终找到的最优位置完全相同。**IG 只有在带剪枝（LB Sifting / Linear Sifting）时才可能影响结果**——因为剪枝会提前终止某个方向，此时先探索哪个方向直接决定能否发现更优位置。

---

## 一句话总结

> **Interaction Graph = 量子电路的 qubit 协作关系图。用它指导 sifting：让合作频繁的 qubit 在 DD 中做邻居。**
