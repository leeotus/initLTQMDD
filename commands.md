# Choi 矩阵 + igGroupSifting 对比实验

## 实验背景

**量子过程层析（QPT）** 通过 Choi 矩阵 `Λ = U ⊗ U*` 完整刻画量子过程。
Choi 矩阵是 2n-qubit 系统（原始电路 n qubit），空间复杂度：

| 表示方式 | 空间 | n=8 示例 |
|----------|------|----------|
| Dense（Qiskit/QuTiP）| `16^n × 16 B` | 281 TB |
| DD（本项目）| O(nodes) 节点数 | ~KB 级 |
| DD + igGroupSifting | O(nodes) 更小 | 视对称性进一步压缩 |

Choi 矩阵具有天然的上下子系统对称性：qubit i 与 qubit i+n 是 U 和 U* 的镜像。
igGroupSifting 通过 InteractionGraph 感知这种对称性，在 Sifting 基础上做对称分组优化。

---

## 对比实验设计

### Exp-1：DD vs Dense 空间对比（motivation）

- **对比对象**：DD 表示 vs Dense 矩阵
- **度量**：DD 节点数 vs `4^(2n)` 元素数
- **目的**：量化 DD 框架本身的空间优势，作为后续实验的 baseline

### Exp-2：三种压缩策略对比（igGroupSifting 的增量贡献）

- **对比对象**：None / Sifting / IGGroupSifting 三种变量序策略
- **度量**：压缩后节点数 / 原始节点数（ig_ratio、sift_ratio，越小越好）
- **假设**：igGroupSifting 对 Clifford 等随机电路 Choi 矩阵更有效；在 Sifting 会破坏结构的场景（如 Grover）igGroupSifting 更稳定

### Exp-3：端到端 QPT 流水线对比（核心对比实验）

- **任务**：Build Choi → 压缩 → partial trace 得到约化密度矩阵
- **对比对象**：三条流水线
  - `None`：直接做 partial trace，无压缩开销
  - `Sifting`：先 Sifting 压缩 Choi，再做 partial trace
  - `IGGroupSifting`：先 igGroupSifting 压缩 Choi，再做 partial trace
- **度量**：
  - `peak_nodes`：流水线中峰值节点数（Choi 压缩后 vs partial trace 结果，取大者）
  - `total_ms`：压缩耗时 + partial trace 耗时
- **意义**：在真实 QPT 任务中，是否应该在 partial trace 之前先压缩？哪种策略的综合代价最小？

### Exp-4：不同门类型横向对比

- **对比对象**：H_layer（高对称）/ CNOT_chain（线性纠缠）/ QFT（规则）/ Grover（不规则）/ Clifford（随机）
- **度量**：各策略在不同门类型上的 peak_nodes 和 total_ms
- **假设**：igGroupSifting 的收益依赖电路的对称性，并非对所有门类型都优于 Sifting

---

## 编译

```bash
cd /home/lijianxian/workshop/initLTQMDD
mkdir -p build_choi_release && cd build_choi_release
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_QFR_TESTS=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -Wno-dev
make benchmark_choi_matrix -j$(nproc)
```

---

## 运行实验

```bash
# Exp-1 & Exp-2 & Exp-4：压缩率对比（compare 模式，默认）
bash test/run_choi_experiments.sh --nmin 2 --nmax 8
# 或直接调用：
./build_choi_release/test/benchmark_choi_matrix --nmin 2 --nmax 8

# Exp-3：端到端 QPT 流水线对比（e2e 模式）
./build_choi_release/test/benchmark_choi_matrix --mode e2e --nmin 2 --nmax 8 \
  > results/choi_matrix/choi_e2e_latest.csv
```

结果文件：
- `results/choi_matrix/choi_results_latest.csv`（compare 模式）
- `results/choi_matrix/choi_e2e_latest.csv`（e2e 模式）

---

## 查看实验结果

### Exp-1：DD 节点数 vs Dense 大小

```bash
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' '$8>0 {ratio=$8/$3;
    printf "gate=%-12s n=%s  dd_nodes=%-6s  dense_elems=%-12s  ratio=%.0fx\n",
    $1,$2,$3,$8,ratio}' | sort -k2,2n -k1,1
```

### Exp-2：igGroupSifting vs Sifting vs None 压缩率

```bash
# 按 ig_ratio 排序
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' 'BEGIN{print "gate,n,nodes_none,nodes_sift,nodes_ig,ig_ratio,sift_ratio"}
    {printf "%s,%s,%s,%s,%s,%s,%s\n",$1,$2,$3,$4,$5,$6,$7}' | \
  sort -t',' -k6 -n | column -t -s','

# Grover 专项：Sifting 反而膨胀（sift_ratio > 1）
awk -F',' 'NR>1 && $1=="Grover" {
  flag=($7+0>1)?"← Sifting膨胀！":"";
  printf "n=%-2s  ig_ratio=%-6s  sift_ratio=%-6s  %s\n",$2,$6,$7,flag}' \
  results/choi_matrix/choi_results_latest.csv

# 胜负统计
awk -F',' 'NR>1 {ig=$6+0; sift=$7+0;
  if(ig<sift) iw++; else if(ig==sift) tie++; else sw++}
  END{print "ig_wins="iw"  tie="tie"  sift_wins="sw}' \
  results/choi_matrix/choi_results_latest.csv
```

### Exp-3：端到端 QPT 流水线对比（核心）

```bash
# 表格形式查看全部
column -t -s',' results/choi_matrix/choi_e2e_latest.csv

# 对每个 (gate, n) 对比三条流水线的 peak_nodes 和 total_ms
awk -F',' 'NR>1 {print $1","$2","$3","$6","$9}' \
  results/choi_matrix/choi_e2e_latest.csv | \
  awk -F',' 'BEGIN{print "gate,n,strategy,peak_nodes,total_ms"} {print}' | \
  column -t -s','

# igGroupSifting vs None：peak_nodes 节省（只看有对称性的场景）
awk -F',' 'NR>1' results/choi_matrix/choi_e2e_latest.csv | \
  awk -F',' '{key=$1","$2; data[key","$3]=$0}
  END{
    for(k in data) {
      split(data[k],a,",");
      gate=a[1]; n=a[2]; strat=a[3];
      none_peak[gate","n]= (strat=="None") ? a[6]+0 : none_peak[gate","n];
      ig_peak[gate","n]  = (strat=="IGGroupSifting") ? a[6]+0 : ig_peak[gate","n];
    }
    print "gate,n,peak_none,peak_ig,saved";
    for(k in none_peak) if(ig_peak[k]!="")
      printf "%s,%s,%s,%.0f%%\n",k,none_peak[k],ig_peak[k],
        (none_peak[k]-ig_peak[k])*100/none_peak[k];
  }' | column -t -s','

# 找 igGroupSifting 在端到端任务中 total_ms 优于 None 的场景
awk -F',' 'NR>1' results/choi_matrix/choi_e2e_latest.csv | sort -t',' -k1,1 -k2,2n -k3,3 | \
  awk -F',' 'BEGIN{print "gate,n,none_total_ms,ig_total_ms,verdict"}
  {if($3=="None") none_ms[$1","$2]=$9;
   if($3=="IGGroupSifting") ig_ms[$1","$2]=$9}
  END{for(k in none_ms) {
    verdict=(ig_ms[k]+0 < none_ms[k]+0)?"ig_faster":"none_faster_or_tie";
    printf "%s,%s,%s,%s\n",k,none_ms[k],ig_ms[k],verdict}}' | \
  sort -t',' -k4,4 -k1,1 | column -t -s','
```

### Exp-4：不同门类型汇总

```bash
# 各门类型的 e2e peak_nodes（三策略对比）
awk -F',' 'NR>1 {print $1","$2","$3","$6}' results/choi_matrix/choi_e2e_latest.csv | \
  awk -F',' 'BEGIN{print "gate,n,strategy,peak_nodes"} {print}' | \
  column -t -s','
```

---

## 核心结论（n=8 实测数据）

| gate | 策略 | peak_nodes | total_ms | 结论 |
|------|------|-----------|---------|------|
| **Grover** | None | 73 | 954 ms | baseline |
| **Grover** | Sifting | **150** | 1345 ms | 峰值膨胀 2×，总时间 +41% ❌ |
| **Grover** | IGGroupSifting | 73 | 1040 ms | 峰值不变，比 Sifting 少 23% ✓ |
| **Clifford** | None | 595 | 505 ms | baseline |
| **Clifford** | Sifting | 403 | 789 ms | 峰值 -32%，但总时间 +56% |
| **Clifford** | IGGroupSifting | 403 | 854 ms | 峰值 -32%，总时间 +69% |
| **QFT** | None | 851 | 485 ms | baseline |
| **QFT** | Sifting | **511** | 1691 ms | 峰值 -40%，但总时间 +249% ❌ |
| **QFT** | IGGroupSifting | 851 | **573 ms** | 峰值不变，总时间仅 +18% ✓ |

**igGroupSifting 的核心价值**：
1. **防御性**：在 Sifting 会破坏结构（Grover）时保持稳定，不引入峰值膨胀
2. **低开销**：即使不产生压缩（QFT），igGroupSifting 的压缩阶段耗时（~99ms）远小于 Sifting（~1108ms）
3. **对称性感知**：对 Clifford 等随机电路能实现与 Sifting 相同的压缩效果（peak -32%）

---

## CSV 列说明

### compare 模式（choi_results_latest.csv）

| 列 | 含义 |
|----|------|
| `gate` | 量子门类型 |
| `n` | 原始电路 qubit 数（Choi 系统为 2n qubit） |
| `choi_nodes_none/sift/ig` | 三种策略下压缩后节点数 |
| `ig_ratio / sift_ratio` | 压缩率（<1 有压缩，>1 反而膨胀） |
| `dense_elems / dense_mb` | Dense 理论大小 |
| `pt_after_none/ig` | partial trace 后节点数 |
| `pt_inflation_*` | partial trace 中间膨胀量 |

### e2e 模式（choi_e2e_latest.csv）

| 列 | 含义 |
|----|------|
| `strategy` | None / Sifting / IGGroupSifting |
| `choi_nodes` | 压缩完成后 Choi DD 节点数 |
| `pt_nodes` | partial trace 结果节点数 |
| `peak_nodes` | max(choi_nodes, pt_nodes)，流水线峰值空间 |
| `compress_ms` | 压缩阶段耗时（None=0） |
| `pt_ms` | partial trace 耗时 |
| `total_ms` | 整个 QPT 流水线总耗时 |
