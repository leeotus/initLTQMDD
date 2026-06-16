#!/usr/bin/env bash
# =============================================================================
# run_dm_experiments.sh  —  密度矩阵噪声仿真对比实验
#
# 用法:
#   bash run_dm_experiments.sh [选项]
#
# 选项:
#   --nmin N            最小 qubit 数 (默认: 4)
#   --nmax N            最大 qubit 数 (默认: 12)
#   --noise-p P         噪声强度 (默认: 0.05)
#   --rounds K          rounds 实验的噪声轮数 (默认: 10)
#   --seeds "s1,s2,..."  Clifford 电路随机种子 (默认: "42,123")
#   --circuits DIR      外部电路目录 (.real/.qasm)，可选
#   --experiments "e1 e2 ..."  运行哪些实验 (默认: "scale rounds compress expval")
#   --output DIR        结果输出目录 (默认: results/density_matrix)
#   --build-type TYPE   Release 或 Debug (默认: Release)
#   -h, --help          显示帮助
#
# 示例:
#   bash run_dm_experiments.sh --nmin 4 --nmax 15 --noise-p 0.1
#   bash run_dm_experiments.sh --nmin 4 --nmax 10 --circuits ~/workshop/circuits
#   bash run_dm_experiments.sh --experiments "scale compress" --nmax 8
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# 默认参数
# ---------------------------------------------------------------------------
NMIN=4
NMAX=12
NOISE_P=0.05
ROUNDS=10
SEEDS="42,123"
CIRCUITS_DIR=""
EXPERIMENTS="scale rounds compress expval"
OUTPUT_DIR="${PROJECT_ROOT}/results/density_matrix"
BUILD_TYPE="Release"

# ---------------------------------------------------------------------------
# 解析参数
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --nmin)        NMIN="$2";        shift 2 ;;
        --nmax)        NMAX="$2";        shift 2 ;;
        --noise-p)     NOISE_P="$2";     shift 2 ;;
        --rounds)      ROUNDS="$2";      shift 2 ;;
        --seeds)       SEEDS="$2";       shift 2 ;;
        --circuits)    CIRCUITS_DIR="$2";shift 2 ;;
        --experiments) EXPERIMENTS="$2"; shift 2 ;;
        --output)      OUTPUT_DIR="$2";  shift 2 ;;
        --build-type)  BUILD_TYPE="$2";  shift 2 ;;
        -h|--help)
            sed -n '2,30p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# 颜色输出
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'
info()    { echo -e "${BLUE}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[ OK ]${NC}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERR ]${NC}  $*"; exit 1; }
section() { echo -e "\n${CYAN}══════════ $* ══════════${NC}"; }

# ---------------------------------------------------------------------------
# Step 1: 确定 build 目录和 binary
# ---------------------------------------------------------------------------
if [[ "${BUILD_TYPE}" == "Release" ]]; then
    BUILD_DIR="${PROJECT_ROOT}/build_dm_release"
else
    BUILD_DIR="${PROJECT_ROOT}/build_dm"
fi
BINARY="${BUILD_DIR}/test/benchmark_density_matrix"

section "检查 ${BUILD_TYPE} Build"
info "Build 目录: ${BUILD_DIR}"

if [ ! -f "${BINARY}" ]; then
    warn "未找到可执行文件，开始编译..."
    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"
    cmake "${PROJECT_ROOT}" \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -Wno-dev \
        > cmake_configure.log 2>&1 \
        || error "cmake 配置失败，查看: ${BUILD_DIR}/cmake_configure.log"
    info "cmake 配置完成，编译 benchmark_density_matrix..."
    make benchmark_density_matrix -j"$(nproc)" \
        > cmake_build.log 2>&1 \
        || error "编译失败，查看: ${BUILD_DIR}/cmake_build.log"
    cd "${PROJECT_ROOT}"
    success "编译完成: ${BINARY}"
else
    # binary 已存在但源码可能更新了，检查时间戳
    SRC="${SCRIPT_DIR}/benchmark_density_matrix.cpp"
    if [ "${SRC}" -nt "${BINARY}" ]; then
        warn "源码比 binary 新，重新编译..."
        cd "${BUILD_DIR}"
        make benchmark_density_matrix -j"$(nproc)" > cmake_build.log 2>&1 \
            || error "编译失败，查看: ${BUILD_DIR}/cmake_build.log"
        cd "${PROJECT_ROOT}"
        success "重编译完成"
    else
        success "Binary 已是最新: ${BINARY}"
    fi
fi

# ---------------------------------------------------------------------------
# Step 2: 显示配置
# ---------------------------------------------------------------------------
section "实验配置"
echo "  Qubit 范围    : ${NMIN} ~ ${NMAX}"
echo "  噪声强度      : p = ${NOISE_P}"
echo "  噪声轮数      : ${ROUNDS}"
echo "  Clifford 种子 : ${SEEDS}"
echo "  运行实验      : ${EXPERIMENTS}"
if [ -n "${CIRCUITS_DIR}" ]; then
    CIRC_COUNT=$(find "${CIRCUITS_DIR}" -maxdepth 1 \( -name "*.real" -o -name "*.qasm" \) 2>/dev/null | wc -l)
    echo "  外部电路目录  : ${CIRCUITS_DIR} (${CIRC_COUNT} 个文件)"
else
    echo "  外部电路      : 仅合成电路 (QFT + RandomClifford)"
fi
echo "  输出目录      : ${OUTPUT_DIR}"
echo ""

# ---------------------------------------------------------------------------
# Step 3: 公共参数构建
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

COMMON_ARGS="--nmin ${NMIN} --nmax ${NMAX} --noise-p ${NOISE_P} --rounds ${ROUNDS} --seeds ${SEEDS}"
if [ -n "${CIRCUITS_DIR}" ]; then
    COMMON_ARGS="${COMMON_ARGS} --circuits ${CIRCUITS_DIR}"
fi

# ---------------------------------------------------------------------------
# 实验执行函数
# ---------------------------------------------------------------------------
run_experiment() {
    local exp_id="$1"
    local mode="$2"
    local header="$3"
    local outfile="${OUTPUT_DIR}/${exp_id}_${TIMESTAMP}.csv"
    local latest="${OUTPUT_DIR}/${exp_id}_latest.csv"

    section "实验: ${exp_id}"
    echo "${header}" > "${outfile}"

    local t_start=$SECONDS
    # shellcheck disable=SC2086
    "${BINARY}" "${mode}" ${COMMON_ARGS} >> "${outfile}" 2>/dev/null
    local t_end=$SECONDS
    local elapsed=$((t_end - t_start))

    cp "${outfile}" "${latest}"
    local rows
    rows=$(tail -n +2 "${outfile}" | wc -l)
    success "完成 (${elapsed}s, ${rows} 行) -> ${outfile}"
    echo ""

    # 终端表格预览（最多20行）
    echo "  预览 (前20行):"
    column -t -s',' "${outfile}" | head -21 | sed 's/^/  /'
    echo ""
}

# ---------------------------------------------------------------------------
# Step 4: 运行选定实验
# ---------------------------------------------------------------------------
for exp in ${EXPERIMENTS}; do
    case "${exp}" in
        scale)
            run_experiment "exp1_scale" "scale" \
                "source,type,n,rho_pure_nodes,rho_noisy_nodes,noise_p,dense_MB,dense_nodes,dense_to_dd_ratio"
            ;;
        rounds)
            run_experiment "exp2_rounds" "rounds" \
                "source,type,n,round,dd_size,purity,trace,dense_nodes"
            ;;
        compress)
            run_experiment "exp3_compress_rho" "compress" \
                "source,type,n,noise_p,size_noisy,size_ig,size_sift,ig_ratio,sift_ratio,ig_ms,sift_ms"
            ;;
        expval)
            run_experiment "exp4_expval" "expval" \
                "source,type,n,noise_p,expval_Z0,analytical,purity,trace"
            ;;
        *)
            warn "未知实验: ${exp}，跳过"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Step 5: 汇总
# ---------------------------------------------------------------------------
section "实验结果汇总"
echo "  时间戳: ${TIMESTAMP}"
echo ""
printf "  %-45s %8s %6s\n" "文件" "行数" "大小"
printf "  %-45s %8s %6s\n" "----" "----" "----"
for f in "${OUTPUT_DIR}"/*_latest.csv; do
    [ -f "$f" ] || continue
    rows=$(tail -n +2 "$f" | wc -l)
    sz=$(du -sh "$f" | cut -f1)
    printf "  %-45s %8s %6s\n" "$(basename $f)" "${rows}" "${sz}"
done
echo ""
success "所有结果保存在: ${OUTPUT_DIR}/"
