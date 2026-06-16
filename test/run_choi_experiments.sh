#!/usr/bin/env bash
# run_choi_experiments.sh
# Choi 矩阵 + igGroupSifting 上层应用对比实验
#
# 用法：
#   bash test/run_choi_experiments.sh [--nmin N] [--nmax N] [--gate GATES] [--memlimit MB]
#                                    [--output DIR] [--filename NAME]
# 示例：
#   bash test/run_choi_experiments.sh --nmin 2 --nmax 10 --gate QFT,Grover
#   bash test/run_choi_experiments.sh --nmin 2 --nmax 8 --memlimit 32000
#   bash test/run_choi_experiments.sh --nmin 2 --nmax 6 --output /tmp/my_results --filename qft_compare.csv

set -euo pipefail

# ============================================================
# 参数解析
# ============================================================
NMIN=2
NMAX=8
GATES=""
MODE="compare"
MEMLIMIT=0
OUT_DIR_ARG=""
FILENAME_ARG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --nmin)      NMIN="$2";        shift 2 ;;
        --nmax)      NMAX="$2";        shift 2 ;;
        --gate)      GATES="$2";       shift 2 ;;
        --mode)      MODE="$2";        shift 2 ;;
        --memlimit)  MEMLIMIT="$2";    shift 2 ;;
        --output)    OUT_DIR_ARG="$2"; shift 2 ;;
        --filename)  FILENAME_ARG="$2"; shift 2 ;;
        -h|--help)
            echo "用法: bash test/run_choi_experiments.sh [选项]"
            echo ""
            echo "选项:"
            echo "  --nmin N         最小 qubit 数 (默认: 2)"
            echo "  --nmax N         最大 qubit 数 (默认: 8)"
            echo "  --gate GATES     门类型过滤，逗号分隔"
            echo "                    可选: H_layer,CNOT_chain,QFT,Grover,Clifford"
            echo "  --mode MODE      compare | e2e (默认: compare)"
            echo "  --memlimit MB    内存限制 (MB)，0=不限制 (默认: 0)"
            echo "  --output DIR     输出目录 (默认: results/choi_matrix)"
            echo "  --filename NAME  输出文件名 (默认: choi_results_<时间戳>.csv)"
            echo "  -h, --help       显示本帮助"
            echo ""
            echo "示例:"
            echo "  bash test/run_choi_experiments.sh --nmin 2 --nmax 10 --gate QFT,Grover"
            echo "  bash test/run_choi_experiments.sh --nmin 2 --nmax 8 --memlimit 32000"
            echo "  bash test/run_choi_experiments.sh --output /tmp/my_results --filename qft_compare.csv"
            exit 0
            ;;
        *) echo "错误: 未知参数 '$1'，使用 --help 查看帮助"; exit 1 ;;
    esac
done

# ============================================================
# 路径
# ============================================================
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build_choi_release"
BINARY="${BUILD_DIR}/test/benchmark_choi_matrix"
SRC="${ROOT}/test/benchmark_choi_matrix.cpp"

if [[ -n "${OUT_DIR_ARG}" ]]; then
    OUT_DIR="${OUT_DIR_ARG}"
else
    OUT_DIR="${ROOT}/results/choi_matrix"
fi
mkdir -p "${OUT_DIR}"

# ============================================================
# 编译（如需要）
# ============================================================
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

# ============================================================
# 组装命令
# ============================================================
CMD="${BINARY} --nmin ${NMIN} --nmax ${NMAX} --mode ${MODE}"

if [[ -n "${GATES}" ]]; then
    CMD="${CMD} --gate ${GATES}"
fi

if [[ "${MEMLIMIT}" -gt 0 ]]; then
    CMD="${CMD} --memlimit ${MEMLIMIT}"
fi

# ============================================================
# 输出文件
# ============================================================
TS=$(date +%Y%m%d_%H%M%S)

if [[ -n "${FILENAME_ARG}" ]]; then
    # 用户指定了文件名
    CSV_RESULT="${OUT_DIR}/${FILENAME_ARG}"
    # 同时生成一个带时间戳的备份
    BASE="${FILENAME_ARG%.csv}"
    CSV_TS="${OUT_DIR}/${BASE}_${TS}.csv"
else
    # 默认：生成时间戳文件 + latest 快捷文件
    CSV_RESULT="${OUT_DIR}/choi_results_${TS}.csv"
fi

# 总是维护一份 choi_results_latest.csv 作为快捷入口
CSV_LATEST="${OUT_DIR}/choi_results_latest.csv"

# ============================================================
# 运行
# ============================================================
echo ""
echo "========================================"
echo " Choi 矩阵 igGroupSifting 对比实验"
echo " nmin=${NMIN}  nmax=${NMAX}  mode=${MODE}"
if [[ -n "${GATES}" ]]; then echo " gate=${GATES}"; fi
if [[ "${MEMLIMIT}" -gt 0 ]]; then echo " memlimit=${MEMLIMIT} MB"; fi
echo " 输出目录: ${OUT_DIR}"
echo " 输出文件: ${CSV_RESULT}"
echo "========================================"
echo ""

if [[ "${MODE}" == "e2e" ]]; then
    if [[ ! -f "${CSV_RESULT}" ]]; then
        "${BINARY}" --nmin 2 --nmax 2 --mode e2e 2>/dev/null | head -1 > "${CSV_RESULT}"
    fi
    eval "${CMD}" >> "${CSV_RESULT}"
else
    eval "${CMD}" > "${CSV_RESULT}"
fi

# 总是更新 latest 快捷文件
cp "${CSV_RESULT}" "${CSV_LATEST}"

# 如果有时间戳备份，也生成
if [[ -n "${CSV_TS:-}" ]]; then
    cp "${CSV_RESULT}" "${CSV_TS}"
    echo "[备份] ${CSV_TS}"
fi

echo "[输出] ${CSV_RESULT}"
echo "[快捷] ${CSV_LATEST}"
echo ""

# ============================================================
# 预览
# ============================================================
if [[ "${MODE}" != "e2e" ]]; then
    echo "=== 结果预览 ==="
    if command -v column &>/dev/null; then
        column -t -s',' "${CSV_RESULT}"
    else
        cat "${CSV_RESULT}"
    fi

    if [[ "$(wc -l < "${CSV_RESULT}")" -gt 1 ]]; then
        echo ""
        echo "=== igGroupSifting vs Sifting vs None（按 ig_ratio 排序）==="
        tail -n +2 "${CSV_RESULT}" | \
            awk -F',' 'BEGIN{print "gate,n,none,sift,ig,ig_ratio,sift_ratio"}
                       {printf "%s,%s,%s,%s,%s,%s,%s\n", $1,$2,$3,$4,$5,$6,$7}' | \
            sort -t',' -k6 -n | column -t -s','

        echo ""
        echo "=== total time (ms) 对比 ==="
        tail -n +2 "${CSV_RESULT}" | \
            awk -F',' 'BEGIN{printf "%-14s %-3s %8s %8s %8s %8s %8s\n","gate","n","build","sift","ig","pt_none","pt_ig"}
                       {printf "%-14s %-3s %8.1f %8.1f %8.1f %8.1f %8.1f\n", $1,$2,$10,$11,$12,$13,$14}' | \
            column -t -s' '
    fi
fi
