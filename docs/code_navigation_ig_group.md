# IG / Group Sifting 代码导读

本文档帮助你快速定位 Interaction Graph、Group Sifting 及相关策略的代码入口和调用链。

---

## 1. 整体架构

```
用户命令
  └─ apps/main.cpp                     # CLI 入口，解析策略名
       └─ QuantumComputation::buildFunctionalityDynamic()
            ├─ 构造阶段：逐门构建 DD，超阈值触发中间 sifting
            └─ 最终阶段：对完整 DD 做一次完整 sifting
                 └─ Package::igGroupSifting() / igLbSifting() / ...
```

---

## 2. 入口文件

| 文件 | 作用 |
|------|------|
| `apps/main.cpp` | CLI 入口，打印用法，调用 buildFunctionalityDynamic |
| `src/QuantumComputation.cpp:1233` | `buildFunctionalityDynamic()` — 动态构造+sifting 主流程 |
| `src/QuantumComputation.cpp:1215` | `parseStrategy()` — 策略名字符串→枚举映射 |

---

## 3. 策略枚举定义

```
extern/dd_package/include/DDpackage.h:156
```

```cpp
enum DynamicReorderingStrategy {
    None, Sifting, LBSifting, TightLBSifting,
    IGSifting, IGLBSifting,
    GroupSifting, IGGroupSifting,
    // ... linear sifting variants ...
};
```

---

## 4. Interaction Graph

### 4.1 数据结构

```
extern/dd_package/include/InteractionGraph.h
```

核心成员：
- `weight[i][j]` — qubit i 和 j 之间的门交互次数
- `degree[i]` — qubit i 的交互度（邻居数量）
- `symmetricGroups` — 对称组检测结果
- `groupId[i]` — qubit i 所属的对称组 ID（-1=单独）

关键方法：
- `addGate(op)` — 每施加一个门时更新图（在构造循环中调用）
- `detectSymmetry()` — 基于 profile hash 检测对称 qubit 组
- `shouldSiftUpFirst()` — 根据 IG 重心决定 sifting 方向

### 4.2 IG 在构造流程中的使用

```cpp
// src/QuantumComputation.cpp:1261
for (auto& op : ops) {
    ig.addGate(op);       // 每个门都更新 IG
    // ... 构建 DD ...
    // 超阈值时用 ig 辅助 sifting
}
```

---

## 5. Sifting 算法实现

### 5.1 基础 Sifting

```
extern/dd_package/src/DDreorder.cpp
```

| 函数 | 行号 | 说明 |
|------|------|------|
| `Package::sifting()` | ~497 | 标准 sifting：最大节点优先选变量，遍历所有位置 |
| `Package::lbSifting()` | 后续 | 带 Lower Bound 剪枝的 sifting |
| `Package::tightLbSifting()` | 后续 | 带 Tight LB 剪枝 |
| `Package::exchangeBaseCase(i, in)` | ~359 | 相邻层交换（不更新 varMap） |
| `Package::exchangeBaseCase2(p, i, in)` | ~398 | 单节点级别的交换操作 |
| `Package::dynamicReorder()` | ~460 | 策略分发入口 |

### 5.2 IG Sifting

```
extern/dd_package/src/DDigSifting.cpp
```

| 函数 | 说明 |
|------|------|
| `Package::igSifting()` | IG 决定变量顺序+方向，无剪枝 |
| `Package::igLbSifting()` | IG 方向 + LB 剪枝（核心改进策略） |
| `selectNextVariable()` | 静态辅助：IG 加权选下一个 sift 变量 |
| `shouldSiftUpFirst()` | 静态辅助：IG 重心决定方向 |

### 5.3 Group Sifting

```
extern/dd_package/src/DDgroupSifting.cpp
```

| 函数 | 说明 |
|------|------|
| `Package::groupSifting()` | 对称检测 + 标准 sifting + 组放置 |
| `Package::igGroupSifting()` | Phase1: igLbSifting → Phase2: 保守组微调 |

---

## 6. Lower Bound 计算

```
extern/dd_package/src/DDlinear.cpp
```

| 函数 | 说明 |
|------|------|
| `computeLowerBoundDown(varMap, i)` | 向下移动时的 LB 估计 |
| `computeLowerBoundUp(varMap, i)` | 向上移动时的 LB 估计 |
| `computeTightLowerBoundDown(...)` | Tight LB（更精确但更慢） |
| `computeTightLowerBoundUp(...)` | Tight LB 向上版本 |
| `countDistinctChildren(level)` | LB 辅助：统计不同子节点数 |

---

## 7. 相邻层交换（Exchange）

这是所有 sifting 的原子操作。有两个重载：

| 签名 | 位置 | 区别 |
|------|------|------|
| `exchangeBaseCase(i, in)` | `DDreorder.cpp:359` | 只交换 DD 结构，不更新 varMap |
| `exchangeBaseCase(i, in, varMap)` | `DDlinear.cpp:353` | 交换 DD + 更新 varMap + 记录 opSeq |

Group Sifting 使用带 varMap 的版本。

---

## 8. Linear Sifting（扩展阅读）

```
extern/dd_package/src/DDlinearv2.cpp  — 主要实现
extern/dd_package/src/DDlinearv3.cpp  — 变体
extern/dd_package/src/DDlinearv4.cpp  — 变体
extern/dd_package/src/DDlinear.cpp    — 基础工具函数
```

Linear Sifting 在普通 exchange 基础上增加 XOR 线性变换，用 `linearInPlace()` 操作。

---

## 9. 调用关系图

```
buildFunctionalityDynamic (src/QuantumComputation.cpp:1233)
│
├─ [构造循环] 逐门 multiply，超阈值触发：
│   ├─ igLbSifting()        ← iggroup/group 策略的中间触发
│   ├─ igSifting()          ← ig 策略
│   ├─ lbSifting()          ← lb 策略
│   └─ sifting()            ← 默认
│
└─ [最终 pass] 构造完成后：
    ├─ igGroupSifting()     ← iggroup 策略
    │   ├─ igLbSifting()         # Phase 1：获取好的基础排序
    │   └─ group adjustment      # Phase 2：对称成员微调
    ├─ groupSifting()       ← group 策略
    │   ├─ sifting 主循环
    │   └─ group placement
    ├─ igLbSifting()        ← iglb 策略
    │   └─ selectNextVariable() + shouldSiftUpFirst() + LB pruning
    └─ sifting()            ← 默认
```

---

## 10. 建议阅读顺序

1. **`InteractionGraph.h`** — 理解 IG 数据结构和 `detectSymmetry()`
2. **`DDreorder.cpp` 的 `sifting()`** — 理解标准 sifting 流程
3. **`DDigSifting.cpp` 的 `igLbSifting()`** — 理解 IG 如何指导方向 + LB 剪枝
4. **`DDgroupSifting.cpp`** — 理解对称组检测如何与 sifting 结合
5. **`DDlinear.cpp` 的 `exchangeBaseCase()`** — 理解原子交换操作
6. **`src/QuantumComputation.cpp:1233`** — 理解整体构造+sifting 流程

---

## 11. 关键注意事项

- `exchangeBaseCase` 会调用 `initComputeTable()` 清空缓存，频繁调用代价高
- QMDD 的 swap 操作在理论上可逆，但实践中因节点共享丢失，"试探+撤销"模式可能导致 DD 膨胀
- `igGroupSifting` 已重构为两阶段设计：先用 igLbSifting 获得好排序，再做保守的局部组调整（MAX_MOVE_DIST=2）
- 对称检测依赖完整 IG，中间构造阶段 IG 不完整时不应做 group placement
