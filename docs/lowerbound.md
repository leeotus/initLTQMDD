Lower bound 计算用的是项目里已有的两个函数。看看具体逻辑：## Lower Bound 原理

DD 变量序为 `0, 1, ..., n-1`（0 最底，n-1 最顶）。变量 `i` 当前在位置 `pos`，`active[varMap[j]]` 表示第 `j` 层的节点数。

### `computeLowerBoundDown` — 变量即将向下移动时的下界

把 DD 节点分三区：

| 区域 | 层范围 | 含义 |
|------|--------|------|
| A | `0 ~ i-1` (下方) | 不受影响的层 |
| 当前层 | `i` | 正在 sift 的变量 |
| B | `i+1 ~ n` (上方) | 变量下移后会被"压过"的层 |

公式：
\[
LB_{down} = \underbrace{\sum_{j<i} \text{active}[j]}_{labA} + \max\left(\text{active}[i],\; 1 + \lfloor labB/2 \rfloor\right)
\]

- **labA**：下方各层节点数之和，这些层不会因为变量 `i` 下移而改变，直接加进来
- **第二项**：变量 `i` 的节点不会因下移而减少到比 `active[i]` 更少；同时，上方 B 区的节点在变量序交换后至少保留 `⌈labB/2⌉` 个（最坏情况下每两个节点合并为一个）— 这给出一个不可能更低的估计

### `computeLowerBoundUp` — 变量即将向上移动时的下界

\[
LB_{up} = \underbrace{\sum_{j \le i-2} \text{active}[j]}_{labA} + 1 + \lfloor\text{active}[i]/2\rfloor + \underbrace{\sum_{j>i} \text{active}[j]}_{labC}
\]

- **labA**：更下方的层（`0` 到 `i-2`），不受变量上移影响
- **1 + active[i]/2**：变量 `i` 上移后，其节点最少能合并到一半（加1是因为根到该层至少有一条路径）
- **labC**：上方各层已经固定，不会减少

### 在 sifting 中如何使用

```
while (pos > 0) {
    lb = computeLowerBoundDown(varMap, pos);
    if (lb > min) break;    // 下界已超过当前最优 → 继续下移不可能更好，剪枝
    exchangeBaseCase(pos, in, varMap);
    ...
}
```

**每次交换前**计算 lower bound。若 LB > 已找到的最小 size (`min`)，说明无论怎么继续移动，DD 总节点数不可能比 `min` 更好 → 提前停止该方向。

### 关键点

这是 BDD 论文（Friedman/Supowit 1990, Jeong/Somenzi 2004 等）中 lower bound sifting 思路的 QMDD 适配。BDD 中每个节点恰好 2 条边，QMDD 中每个节点有 4 条边（RADIX=2, NEDGE=4），所以 `/2` 的系数是一个偏保守的估计——实际 QMDD 中合并概率可能更低，导致 lower bound 偏松。这也解释了为什么 LB Sifting 结果（956）略差于完整 Sifting（932）：剪枝跳过了一些本可以找到更优位置的交换。