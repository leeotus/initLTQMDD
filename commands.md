# Choi 矩阵 + igGroupSifting 对比实验

## 实验背景

**量子过程层析（QPT）** 通过 Choi 矩阵 `Λ = U ⊗ U*` 完整刻画量子过程。
Choi 矩阵是 2n-qubit 系统（原始电路 n qubit），空间复杂度：

| 表示方式 | 空间 | n=8 |
|----------|------|-----|
| Dense（Qiskit/QuTiP）| `16^n × 16 B` | **281 TB** |
| DD（本项目，无压缩）| O(nodes) | ~1 KB |
| DD + igGroupSifting | O(nodes) 更小 | 更小 |

Choi 矩阵具有天然的上下子系统对称性：qubit i 与 qubit i+n 是 U 和 U* 的镜像。
igGroupSifting 通过 InteractionGraph 感知这种对称性，在 Sifting 基础上做对称分组优化。

---

## 对比实验设计

### Exp-1：DD 节点数 vs Dense 理论大小

**对比对象**：DD 表示（本项目）vs Dense 矩阵（Qiskit/QuTiP）  
**对比维度**：存储节点数 vs 理论元素数  
**假设**：DD 节点数随 n 线性增长，Dense 指数爆炸

```
DD nodes ~ O(n)   vs   Dense elems = 4^(2n) = 16^n
```

### Exp-2：igGroupSifting vs Sifting vs 无压缩

**对比对象**：三种变量序策略  
**对比维度**：压缩后节点数 / 原始节点数（压缩率，越小越好）  
**假设**：igGroupSifting 利用 Choi 矩阵对称性，比 Sifting 压缩率更高或更稳定

| 策略 | 原理 | 预期优势 |
|------|------|---------|
| None | 默认序 | baseline |
| Sifting | 贪心逐变量交换 | 对 QFT 等规则门有效 |
| igGroupSifting | 对称组检测 + IG 引导 | 对 Clifford 等随机门更稳定，不破坏已有结构 |

### Exp-3：压缩前 vs 压缩后做 partial trace

**对比对象**：在 igGroupSifting 压缩后 vs 直接做 partial trace  
**对比维度**：partial trace 后的节点数（越小说明压缩对下游操作越有益）  
**假设**：压缩后的 Choi DD 在做 partial trace（求约化密度矩阵）时中间节点更少

### Exp-4：不同量子门类型的对比

**对比对象**：H层（高对称）/ CNOT链（线性纠缠）/ QFT（规则）/ Grover（不规则）/ Clifford（随机）  
**对比维度**：各策略压缩率在不同门类型上的差异  
**假设**：igGroupSifting 对高对称性门（H层、Clifford）效果更好；Sifting 对规则门（QFT）更好

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
# 完整运行（n=2~8，自动保存 CSV + 打印摘要）
bash test/run_choi_experiments.sh

# 自定义 qubit 范围
bash test/run_choi_experiments.sh --nmin 2 --nmax 10

# 直接调用，输出到 stdout
./build_choi_release/test/benchmark_choi_matrix --nmin 2 --nmax 8
```

结果保存在 `results/choi_matrix/choi_results_latest.csv`。

---

## 查看实验对比结果

### Exp-1：DD 节点数 vs Dense 大小

```bash
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' '$8>0 {ratio=$8/$3;
             printf "gate=%-12s n=%s  dd_nodes=%-6s  dense_elems=%-20s  ratio=%.0fx\n",
             $1,$2,$3,$8,ratio}' | sort -k2,2n -k1,1
```

预期结果（n=8 Clifford）：
```
gate=Clifford     n=8  dd_nodes=595    dense_elems=...  ratio=...x
```

### Exp-2：igGroupSifting vs Sifting vs None

```bash
# 压缩率对比（按 ig_ratio 从小到大排序）
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' 'BEGIN{print "gate,n,nodes_none,nodes_sift,nodes_ig,ig_ratio,sift_ratio"}
             {printf "%s,%s,%s,%s,%s,%s,%s\n",$1,$2,$3,$4,$5,$6,$7}' | \
  sort -t',' -k6 -n | column -t -s','

# igGroupSifting 有压缩效果（ig_ratio < 1）的行
awk -F',' 'NR>1 && $6+0 < 1.0 {
  printf "gate=%-12s n=%-2s  ig_ratio=%s  sift_ratio=%s\n",$1,$2,$6,$7}' \
  results/choi_matrix/choi_results_latest.csv

# Grover 专项：igGroupSifting 避免 Sifting 破坏（sift_ratio > 1 = 膨胀）
awk -F',' 'NR>1 && $1=="Grover" {
  printf "n=%-2s  ig_ratio=%s  sift_ratio=%s  %s\n",
  $2,$6,$7,($7+0>1?"← Sifting膨胀！igGroupSifting更稳定":"")}' \
  results/choi_matrix/choi_results_latest.csv

# 胜负统计
awk -F',' 'NR>1 {ig=$6+0; sift=$7+0;
  if(ig<sift) ig_win++; else if(ig==sift) tie++; else sift_win++}
  END{print "ig_wins="ig_win"  tie="tie"  sift_wins="sift_win}' \
  results/choi_matrix/choi_results_latest.csv
```

关键发现：
- **Clifford n=6~8**：igGroupSifting 压缩率 0.68~0.70
- **Grover n=8**：Sifting `ratio=2.05`（反而膨胀），igGroupSifting `ratio=1.0`（保持稳定）
- **QFT**：Sifting 始终有效（0.60~0.76），igGroupSifting 无效（1.0）—— 各有专长

### Exp-3：压缩后 partial trace 的节点膨胀

```bash
# partial trace 节点对比：ig压缩版 vs 无压缩版
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' 'BEGIN{print "gate,n,pt_none,pt_ig,saved,saved%"}
             {if($14>0) pct=int(($14-$17)*100/$14); else pct=0;
              printf "%s,%s,%s,%s,%s,%s%%\n",$1,$2,$14,$17,($14-$17),pct}' | \
  column -t -s','

# 重点：partial trace 节省 > 0 的行（压缩确实有益于下游）
awk -F',' 'NR>1 && ($14-$17)>0 {
  printf "gate=%-12s n=%-2s  pt_none=%-5s  pt_ig=%-5s  saved=%s\n",
  $1,$2,$14,$17,($14-$17)}' results/choi_matrix/choi_results_latest.csv
```

关键发现：
- **Clifford n=6**：无压缩 pt=1，ig压缩后 pt=74（结构不同，ig 保持更多信息）
- **Clifford n=7**：无压缩 pt=170，ig压缩后 pt=138（**减少 19%**）
- **Clifford n=8**：无压缩 pt=298，ig压缩后 pt=202（**减少 32%**）
- **QFT n=8**：无压缩 pt=426，ig压缩后 pt=426（持平，ig对QFT无效，一致）

### Exp-4：不同门类型对比

```bash
# 各门类型的平均 ig_ratio
awk -F',' 'NR>1 {sum[$1]+=$6+0; cnt[$1]++}
  END{for(g in sum) printf "gate=%-12s  avg_ig_ratio=%.4f\n",g,sum[g]/cnt[g]}' \
  results/choi_matrix/choi_results_latest.csv | sort -k2

# 各门类型的平均 sift_ratio
awk -F',' 'NR>1 {sum[$1]+=$7+0; cnt[$1]++}
  END{for(g in sum) printf "gate=%-12s  avg_sift_ratio=%.4f\n",g,sum[g]/cnt[g]}' \
  results/choi_matrix/choi_results_latest.csv | sort -k2

# 完整对比矩阵（每种门 × 每种策略）
tail -n +2 results/choi_matrix/choi_results_latest.csv | \
  awk -F',' '{print $1","$2","$6","$7}' | \
  awk -F',' 'BEGIN{print "gate,n,ig_ratio,sift_ratio"} {print}' | \
  column -t -s','
```

---

## CSV 列说明

| 列 | 含义 |
|----|------|
| `gate` | 量子门（H_layer/CNOT_chain/QFT/Grover/Clifford） |
| `n` | 原始电路 qubit 数（Choi 系统为 2n qubit） |
| `choi_nodes_none/sift/ig` | 三种策略下 Choi DD 节点数 |
| `ig_ratio` | igGroupSifting 压缩率（节点数之比，<1 有压缩） |
| `sift_ratio` | Sifting 压缩率（>1 表示反而膨胀） |
| `dense_elems` | Dense 所需元素数（4^(2n)） |
| `dense_mb` | Dense 所需内存（MB） |
| `sift_ms / ig_ms` | 各策略耗时（ms） |
| `pt_before_none/ig` | partial trace 前节点数 |
| `pt_after_none/ig` | partial trace 后节点数 |
| `pt_inflation_none/ig` | partial trace 中间膨胀量（负值=有压缩） |
