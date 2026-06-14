# QMDD+QEC vs MQT QECC：详细对比分析

## 1. 概览

| 维度 | **MQT QECC** | **本工具 QMDD+QEC** |
|------|-------------|---------------------|
| **全称** | Munich Quantum Toolkit — Quantum Error Correcting Codes | QMDD-based QEC Simulation Framework |
| **开发者** | TUM (Wille 组) | 本项目 (initLTQMDD + ddsim) |
| **GitHub** | https://github.com/munich-quantum-toolkit/qecc | 本项目 qec/ 模块 |
| **语言** | **纯 Python** | **C++17** |
| **安装** | `pip install mqt.qecc` | `cmake -DBUILD_QEC=ON && make` |
| **许可证** | MIT | MIT |
| **能否被 C++ 调用** | ❌ 不能直接调用（纯 Python 包） | ✅ 原生 C++ 库 |
| **依赖** | Stim, Z3, Qiskit, numpy, scipy, numba, ldpc | dd_package (JKQ/IIC) |
| **定位** | **QEC 理论/设计工具** | **QEC 高效模拟引擎** |
| **工作层** | 逻辑层（码设计、解码、FT 线路合成） | 物理层（编码后电路的 DD 紧凑表示与模拟） |

---

## 2. 功能对比

### 2.1 MQT QECC 的核心功能

| 模块 | 功能 | 方法 |
|------|------|------|
| **`cc_decoder/`** | Color code 解码 | LightOut MaxSAT (Z3 SMT) |
| **`analog_information_decoding/`** | Bosonic QLDPC 模拟解码 | Belief Propagation + OSD with analog info |
| **`circuit_synthesis/`** | Fault-tolerant 状态准备线路合成 | SAT-based / 确定性合成 |
| **`code_switching/`** | 代码切换优化 | 搜索最小切换操作放置 |
| **`cococo/`** | Color code lattice surgery 编译 | CNOT+T 线路 lattice surgery (static/movable qubits) |
| **`codes/`** | 码库 | CSS Code, Surface Code, Color Code, Steane(存为 .npy), Shor, Hamming 等 |

### 2.2 本 QMDD+QEC 的核心功能

| 模块 | 功能 | 方法 |
|------|------|------|
| **`QECCode`** | 编码抽象 + 展开为物理电路 | 虚基类，支持多种编码扩展 |
| **`SteaneCode`** | Steane [[7,1,3]] 完整实现 | 编码/解码/纠错生成 |
| **`QECSimulator`** | 噪声环境下模拟 QEC 电路 | ddsim 噪声引擎 + IG/GroupSifting |
| **`Experiment 1-5`** | 编码验证 / DD 压缩 / 阈值扫描 / 可扩展性 / IG 对称性 | 内置 5 个基准实验 |

---

## 3. 核心差异：工作层次不同

```
抽象层次:
                  ┌──────────────────────────────────────┐
                  │           MQT QECC                   │
   逻辑层         │   · 码设计 (CSS/Surface/Color)       │
 (code design)    │   · FT 状态准备线路合成              │
                  │   · 码距/解码/编码方案               │
                  │   · Lattice surgery 编译             │
                  └───────────┬──────────────────────────┘
                              │ 输出的逻辑线路
                              ▼
                  ┌──────────────────────────────────────┐
   物理层         │       本工具 QMDD+QEC                │
 (simulation)     │   · DD 紧凑表示 QEC 电路             │
                  │   · IG 对称组检测与 Group Sifting    │
                  │   · 高效率噪声模拟                    │
                  │   · DD 可扩展性分析                  │
                  └──────────────────────────────────────┘
```

**互补关系**：
- MQT QECC 设计 QEC 方案 → 输出逻辑线路 → 本工具展开为物理电路并高效模拟
- 两者的输入/输出天然可衔接

---

## 4. 技术路线对比

### 4.1 MQT QECC 的技术路线

```
StabilizerCode (Python)
  ├── 依赖 Stim (stabilizer 快速模拟)
  ├── 依赖 Z3 SMT Solver (MaxSAT 解码)
  ├── 依赖 numpy/scipy (数值计算)
  ├── 依赖 Qiskit (电路构建)
  ├── 依赖 numba (JIT 加速)
  └── 依赖 ldpc (LDPC 码操作)

核心算法:
  - 解码: LightsOut → MaxSAT formula → Z3 Optimize
  - 合成: SAT-based FT state preparation
  - 编译: Lattice surgery in color codes
```

### 4.2 本 QMDD+QEC 的技术路线

```
QECCode → SteaneCode (C++)
  ├── 依赖 JKQ DD Package (QMDD 表示)
  ├── 依赖 InteractionGraph (IG 对称检测)
  ├── 依赖 DDreorder (Group Sifting)
  └── 依赖 ddsim (噪声模拟引擎)

核心算法:
  - 模拟: QMDD 紧凑表示 + DD 乘法
  - 优化: IG 对称组检测 → Group Sifting 压缩 DD
  - 噪声: ddsim 内置 depolarization / APD / phase flip
```

---

## 5. "都做 QEC 模拟"的误解消除

| 看似相同的表述 | MQT QECC 实际做的是 | 本工具实际做的是 |
|-------------|-------------------|-----------------|
| "模拟 QEC" | 模拟**解码器**的纠错能力 | 模拟**编码后的物理电路**的量子态演化 |
| "支持 Steane code" | 存储 stabilizer 矩阵为 `.npy` 文件 | 生成 Steane 码的**完整门序列**并用 DD 紧凑表示 |
| "噪声模拟" | 通过 Stim 的 Pauli frame 注入 syndrome 噪声 | 在**每个物理门**后注入噪声算子，精确影响量子态 |
| "可扩展性" | 更大码距/更多 qubit 的码 | DD size 随 qubit 数和 QEC 轮数的增长趋势 |

---

## 6. 如何在 C++ 项目中使用 MQT QECC？

### 答案：不是直接的

```
MQT QECC 是纯 Python 包，不能在 C++ 中直接 #include。
但可以通过以下方式间接调用:
```

**方式 1：Python subprocess 调用**（推荐用于对比实验）

```cpp
// qec_benchmark 中调用 MQT QECC 获取 Steane stabilizer 矩阵
#include <cstdlib>
#include <fstream>

std::string runPythonScript(const std::string& script) {
    std::string cmd = "python3 -c '" + script + "' > /tmp/mqt_output.txt";
    system(cmd.c_str());
    std::ifstream f("/tmp/mqt_output.txt");
    return std::string(std::istreambuf_iterator<char>(f),
                       std::istreambuf_iterator<char>());
}

// 获取 Steane 码的 HX/HZ 矩阵
void getSteaneFromMQT() {
    std::string script = 
        "import numpy as np; "
        "hx=np.load('path/to/steane/hx.npy'); "
        "print(hx.tolist())";
    std::string result = runPythonScript(script);
    // 解析 JSON 格式的矩阵...
}
```

**方式 2：Python 绑定**（pybind11，本项目的 `extern/pybind11` 已就绪）

```cpp
// 可以通过 pybind11 嵌入 Python 解释器
#include <pybind11/embed.h>
namespace py = pybind11;

void callMQT() {
    py::scoped_interpreter guard{};
    py::module_ mqt = py::module_::import("mqt.qecc.codes");
    auto steane_hx = mqt.attr("np").attr("load")("steane/hx.npy");
    // ...
}
```

**方式 3：离线导出 → 在 C++ 中读取**（最简单）

```bash
# Step 1: 用 Python 导出 MQT QECC 的 stabilizer 数据
python3 -c "
from mqt.qecc.codes import CSSCode
import numpy as np

# Steane code from MQT
hx = np.load('path/to/steane/hx.npy')
hz = np.load('path/to/steane/hz.npy')

# 写入为 C++ 可读的文本格式
np.savetxt('/tmp/steane_hx.txt', hx, fmt='%d')
np.savetxt('/tmp/steane_hz.txt', hz, fmt='%d')
print('Exported to /tmp/steane_hx.txt and /tmp/steane_hz.txt')
"
```

```cpp
// Step 2: 在 C++ 中读取
#include <fstream>
std::vector<std::vector<int>> loadStabMatrix(const std::string& path) {
    std::ifstream f(path);
    std::vector<std::vector<int>> mat;
    int val;
    // parse line by line...
    return mat;
}

auto HX = loadStabMatrix("/tmp/steane_hx.txt");
auto HZ = loadStabMatrix("/tmp/steane_hz.txt");
// 用 HX/HZ 构造 StabilizerCode...
```

---

## 7. 如何做两者的对比实验？

### 7.1 可比维度分析

由于两工具工作在不同层次，**直接对比需聚焦在重叠部分**：

| 对比维度 | 可比性 | 说明 |
|----------|--------|------|
| **Clifford QEC 模拟速度** | ✅ 可比 | 两者都支持 Clifford 电路 (MQT 用 Stim，我们用 DD) |
| **Logical Error Rate** | ✅ 可比 | 同一噪声模型下的 LER |
| **可扩展性** | ⚠️ 需统一定义 | MQT 看 code distance scaling；我们看 DD size scaling |
| **Non-Clifford 门模拟** | ❌ 不可比 | MQT 的 Stim 后端不支持 T 门 |
| **DD 压缩率** | ❌ 不可比 | MQT 无 DD 表示 |
| **解码器性能** | ❌ 不可比 | 我们不实现解码器 |
| **IG 对称性** | ❌ 不可比 | MQT 无 IG 分析 |

### 7.2 可行的对比方案

**方案 A：Clifford QEC 模拟速度对比 (最直接)**

```
公共输入: Steane [[7,1,3]] 编码后的物理电路 (13 qubit)
噪声模型: Depolarization p=0.001
指标: wall-clock 时间 (10 轮 QEC)

MQT QECC 侧:
  from mqt.qecc.codes.steane import steane_code
  import stim
  # 构建 circuit → stim.compile_sampler() → 计时

本工具侧:
  ./build_qec/qec/qec_benchmark
  # 自动输出 LER + DD size + 时间
```

**方案 B：LER 准确性交叉验证**

```
目的: 验证两个工具在相同噪声下的 LER 是否一致
方法:
  1. 用本工具在 depolarization p=0.001 下运行 10000 shots
  2. 用 MQT QECC + Stim 在相同条件下运行 10000 shots
  3. 比较 LER 置信区间是否重叠

预期: LER 在统计误差范围内一致
```

**方案 C：互补验证 (展示两工具的协同价值)**

```
目的: 展示两工具如何互补而非竞争
方法:
  1. MQT QECC: 导出 stabilizer 矩阵
  2. 本工具: 从 stabilizer 生成物理电路并 DD 模拟
  3. MQT QECC: 对模拟结果用 MaxSAT 解码
  4. 完整流程的端到端时间 vs 单独用 Stim 或单独用 Aer

这展示了 1+1>2 的协同效应
```

### 7.3 推荐对比脚本

```python
#!/usr/bin/env python3
"""
compare_mqt_vs_ours.py
对比 MQT QECC (Stim backend) vs Our QMDD+QEC
"""

import subprocess, time, json, csv, os
import numpy as np

# ===== MQT QECC side =====
try:
    from mqt.qecc.codes.stabilizer_code import StabilizerCode
    import stim
    HAS_MQT = True
except ImportError:
    HAS_MQT = False
    print("MQT QECC not installed. Run: pip install mqt.qecc")

# ===== Our tool side =====
OUR_BIN = "./build_qec/qec/qec_benchmark"

def run_mqt(num_rounds, noise_p, shots=1000):
    """Run MQT QECC + Stim on Steane code"""
    if not HAS_MQT:
        return {"status": "SKIP"}
    
    # Build Steane code circuit
    # MQT stores stabilizer as numpy files
    import importlib.resources
    hx = np.load("path/to/steane/hx.npy")
    hz = np.load("path/to/steane/hz.npy")
    
    # Build stim circuit (Clifford gates only for fair comparison)
    # ... circuit construction ...
    
    t0 = time.time()
    # sampler = circuit.compile_sampler()
    # results = sampler.sample(shots)
    elapsed = time.time() - t0
    
    return {"time": elapsed, "ler": 0.0}

def run_ours(num_rounds, noise_p, shots=1000):
    """Run our QMDD+QEC tool"""
    t0 = time.time()
    result = subprocess.run(
        [OUR_BIN], 
        capture_output=True, text=True,
        timeout=300
    )
    elapsed = time.time() - t0
    # Parse output for LER, DD size
    return {"time": elapsed, "ler": 0.0}

def main():
    configs = [
        {"rounds": 1, "noise": 0.0},
        {"rounds": 5, "noise": 0.001},
        {"rounds": 10, "noise": 0.01},
    ]
    
    results = []
    for cfg in configs:
        r_mqt = run_mqt(cfg["rounds"], cfg["noise"])
        r_ours = run_ours(cfg["rounds"], cfg["noise"])
        
        results.append({
            "rounds": cfg["rounds"],
            "noise": cfg["noise"],
            "mqt_time": r_mqt.get("time", "N/A"),
            "ours_time": r_ours.get("time", "N/A"),
            "mqt_ler": r_mqt.get("ler", "N/A"),
            "ours_ler": r_ours.get("ler", "N/A"),
        })
    
    # Write CSV
    with open("comparison_mqt_vs_ours.csv", "w") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    
    print("Results saved to comparison_mqt_vs_ours.csv")

if __name__ == "__main__":
    main()
```

---

## 8. 各自的独特优势 (USP)

### MQT QECC 的 USP

- ✅ **解码器设计**：MaxSAT 解码 / BP+OSD 模拟解码
- ✅ **FT 线路合成**：自动生成容错状态准备线路
- ✅ **Lattice surgery**：color code 专用编译器
- ✅ **Python 生态**：`pip install` 开箱即用
- ✅ **论文发表**：PRX Quantum ×2 + ASP-DAC + VLSID

### 本 QMDD+QEC 的 USP

- ✅ **DD 压缩**：QEC 电路 DD size << state vector size
- ✅ **IG 对称性**：自动检测 QEC 电路 qubit 对称结构
- ✅ **Group Sifting**：利用对称性加速 DD 变量序优化
- ✅ **全门集模拟**：支持 Clifford + non-Clifford (Stim 无法)
- ✅ **C++ 性能**：原生高性能，可直接嵌入其他 C++ 项目

---

## 9. 所属生态系统

| 方面 | MQT QECC | 本工具 QMDD+QEC |
|------|----------|-----------------|
| **所属项目** | Munich Quantum Toolkit (MQT) | JKQ QFR / initLTQMDD |
| **生态系统** | MQT 系列 (qmap, qecc, qcec, qcompiler) | DD-based 量子计算工具链 |
| **核心论文** | Berent et al. Quantum 2024; Schmid et al. PRX Quantum 2025 | Zulehner & Wille TCAD 2018; Grurl et al. TCAD 2023 |
| **上游依赖** | Stim, Z3, numpy, scipy, numba, ldpc, Qiskit | dd_package (JKQ/IIC) |
| **Python API** | ✅ 原生 | ❌ (C++ only, pybind11 planned) |
| **C++ 可用性** | ❌ | ✅ |

---

## 10. 引用

```bibtex
@article{berent2024decoding,
  title   = {Decoding quantum color codes with {MaxSAT}},
  author  = {Berent, Lucas and Burgholzer, Lukas and Derks, Peter-Jan and Eisert, Jens and Wille, Robert},
  journal = {Quantum}, volume = {8}, pages = {1506}, year = {2024},
}

@article{schmid2025deterministic,
  title   = {Deterministic Fault-Tolerant State Preparation for Near-Term Quantum Error Correction},
  author  = {Schmid, Ludwig and Peham, Tom and Berent, Lucas and M{\"u}ller, Markus and Wille, Robert},
  journal = {PRX Quantum}, volume = {6}, pages = {020330}, year = {2025},
}

@inproceedings{grurl2023automatic,
  title     = {Automatic Implementation and Evaluation of Error-Correcting Codes for Quantum Computing},
  author    = {Grurl, Thomas and Pichler, Christoph and Fuss, J{\"u}rgen and Wille, Robert},
  booktitle = {VLSID}, year = {2023},
}
```

---

## 11. 总结

| 问题 | 答案 |
|------|------|
| **MQT QECC 能在 C++ 中使用吗？** | **不能直接使用。** 它是纯 Python 包。需通过 Python subprocess / pybind11 嵌入 / 离线导出数据间接调用 |
| **两工具是竞争关系吗？** | **不是。** 工作在不同层次：MQT 做逻辑层设计，我们做物理层模拟 |
| **如何做对比？** | 在 **Clifford QEC 模拟速度**和 **LER 准确性** 两个维度可比；其他维度（DD 压缩、Non-Clifford、IG 对称）是我们的独有优势 |
| **推荐的论文策略？** | 不是 "vs"，而是 "complementary"——展示 MQT 设计 + 我们模拟的协同 workflow |