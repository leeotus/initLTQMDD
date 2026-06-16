#!/usr/bin/env bash
# run_choi_experiments.sh
# Choi 矩阵 + igGroupSifting 上层应用对比实验
# 实验内容：
#   1. 多种量子门 Choi 矩阵 DD 节点数 vs Dense 理论大小
#   2. igGroupSifting / Sifting / None 三组压缩率对比
#   3. 压缩后做 partial trace 的中间节点膨胀量对比
# 用法：
#   bash test/run_choi_experiments.sh [--nmin N] [--nmax N]
# 示例：
#   bash test/run_choi_experiments.sh --nmin 2 --nmax 8

set -euo pipefail

# ---------- 参数解析 ----------
NMIN=2
NMAX=8
while [[ $# -gt 0 ]]; do
    case "$1" in
        --nmin)  NMIN="$2";  shift 2 ;;
        --nmax)  NMAX="$2";  shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

# ---------- 路径 ----------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build_choi_release"
BINARY="${BUILD_DIR}/test/benchmark_choi_matrix"
OUT_DIR="${ROOT}/results/choi_matrix"
SRC="${ROOT}/test/benchmark_choi_matrix.cpp"

mkdir -p "${OUT_DIR}"

# ---------- 编译 ----------
need_build=0
if [[ ! -f "${BINARY}" ]]; then
    echo "[build] benchmark_choi_matrix 不存在，开始编译..."
    need_build=1
elif [[ "${SRC}" -nt "${BINARY}" ]]; then
    echo "[build] 源文件更新，重新编译..."
    need_build=1
fi

if [[ $need_build -eq 1 ]]; then
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"
    cmake "${ROOT}" -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_QFR_TESTS=ON \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -Wno-dev \
        > cmake_configure.log 2>&1
    make benchmark_choi_matrix -j"$(nproc)" > cmake_build.log 2>&1
    echo "[build] 编译完成: ${BINARY}"
    cd "${ROOT}"
fi

# ---------- 时间戳 ----------
TS=$(date +%Y%m%d_%H%M%S)
CSV_LATEST="${OUT_DIR}/choi_results_latest.csv"
CSV_TS="${OUT_DIR}/choi_results_${TS}.csv"

# ---------- 运行实验 ----------
echo ""
echo "========================================"
echo " Choi 矩阵 igGroupSifting 对比实验"
echo " nmin=${NMIN}  nmax=${NMAX}"
echo "========================================"
echo ""

"${BINARY}" --nmin "${NMIN}" --nmax "${NMAX}" | tee "${CSV_LATEST}"
cp "${CSV_LATEST}" "${CSV_TS}"

echo ""
echo "[输出] ${CSV_LATEST}"
echo "[备份] ${CSV_TS}"

# ---------- 漂亮预览 ----------
echo ""
echo "=== 结果预览 ==="
if command -v column &>/dev/null; then
    column -t -s',' "${CSV_LATEST}"
else
    cat "${CSV_LATEST}"
fi

# ---------- 关键对比摘要 ----------
echo ""
echo "=== igGroupSifting vs Sifting vs None（按 ig_ratio 排序）==="
tail -n +2 "${CSV_LATEST}" | \
    awk -F',' 'BEGIN{print "gate,n,none,sift,ig,ig_ratio,sift_ratio"}
               {printf "%s,%s,%s,%s,%s,%s,%s\n", $1,$2,$3,$4,$5,$6,$7}' | \
    sort -t',' -k6 -n | column -t -s','

echo ""
echo "=== partial trace 膨胀：ig压缩 vs 无压缩（膨胀越小越好）==="
tail -n +2 "${CSV_LATEST}" | \
    awk -F',' 'BEGIN{print "gate,n,pt_inflation_none,pt_inflation_ig,减少"}
               {diff=$16-$13; printf "%s,%s,%s,%s,%s\n", $1,$2,$13,$16,diff}' | \
    column -t -s','

echo ""
echo "=== DD 节点数 vs Dense 元素数（压缩倍数）==="
tail -n +2 "${CSV_LATEST}" | \
    awk -F',' '$8>0 {ratio=$8/$3; printf "gate=%-12s n=%s  dd_nodes=%s  dense_elems=%s  ratio=%.0fx\n",
               $1,$2,$3,$8,ratio}' | sort -k2,2n -k1,1
