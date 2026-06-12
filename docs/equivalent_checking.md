有几个主要的可对比项目：

## 直接可对比的项目

### 1. QCEC（JKQ 团队自己的）⭐ 最相关

- **仓库**: https://github.com/cda-tum/qcec
- **关系**: 与你的项目同源（都基于 JKQ DD package）
- **方法**: DD-based equivalence checking + 多种优化（alternating scheme, gate-cost heuristic）
- **对比点**: 你的 igGroupSifting 加速 vs 他们的 alternating construction scheme
- **论文**: Burgholzer & Wille, "Advanced Equivalence Checking for Quantum Circuits", IEEE TCAD 2021

### 2. QMDD-based Verification（Niemann/Wille 早期工作）

- **论文**: "QMDDs: Efficient Quantum Function Representation and Manipulation", IEEE TCAD 2016
- **方法**: 原始 QMDD 等价性验证，无 sifting 优化
- **对比点**: 你的 IG-guided sifting 对比他们的无优化构建

### 3. ZX-Calculus Based（PyZX）

- **仓库**: https://github.com/Quantomatic/pyzx
- **方法**: 用 ZX-calculus 图重写规则化简，完全不同的理论框架
- **对比点**: 图重写 vs DD 表示，不同电路类型上各有优劣
- **论文**: Kissinger & van de Wetering, "Reducing T-count with the ZX-calculus", 2020

### 4. Quartz（UC Berkeley）

- **仓库**: https://github.com/quantum-compiler/quartz
- **方法**: 基于 ECC（等价电路类）的验证
- **对比点**: 搜索式验证 vs DD 直接计算

## 最推荐对比

| 项目 | 为什么适合对比 | 难度 |
|------|--------------|------|
| **QCEC** | 同框架(DD)、同输入格式(.qasm/.real)、直接跑 benchmark | ⭐ 低 |
| **PyZX** | 完全不同方法论，互补对比有学术价值 | 中 |
| **无sifting基线** | 你自己项目 `none` 策略就是基线 | 零 |

## 实操建议

**QCEC 最容易对比**——它也用 JKQ DD package，输入格式兼容，直接装：

```bash
pip install mqt.qcec
```

Python API：
```python
from mqt import qcec
result = qcec.verify(circ1, circ2)
print(result.equivalence)  # "equivalent" / "not_equivalent"
print(result.check_time)   # 直接拿时间对比
```

你的切入角度：**QCEC 用的是 alternating DD construction（交替从两电路取门），你的是 igGroupSifting + 动态重排序**。在特定电路（高对称性、多 qubit 交互）上你的方法可能更优。