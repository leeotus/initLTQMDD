#!/bin/bash
# 量子电路等价性验证 - 多策略对比实验
# Usage: bash equiv_benchmark.sh [strategy1,strategy2,...]
# 支持逗号分隔多个策略同时对比, 例如: bash equiv_benchmark.sh iggroup,none,sifting
# strategy: none, sifting, lb, tightlb, ig, iglb, group, iggroup (默认 iggroup,none,sifting)

STRATEGIES_INPUT="${1:-iggroup,none,sifting}"
IFS=',' read -ra STRATEGIES <<< "$STRATEGIES_INPUT"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RELEASE_DIR="${PROJECT_DIR}/release"
CIRCUITS_DIR="/home/lijianxian/workshop/circuits"
EQUIV_CHECK="${RELEASE_DIR}/equiv_check"
RESULTS_FILE="${SCRIPT_DIR}/equiv_benchmark_results.csv"

# Build if needed
if [ ! -f "$EQUIV_CHECK" ]; then
    echo "Building equiv_check..."
    cd "$RELEASE_DIR" && make -j$(nproc) equiv_check > /dev/null 2>&1
fi

NUM_STRATS=${#STRATEGIES[@]}

echo "=============================================="
echo " 量子电路等价性验证 - 多策略对比"
echo " 策略: ${STRATEGIES_INPUT} (共${NUM_STRATS}种)"
echo "=============================================="
echo ""

# Build CSV header dynamically
CSV_HEADER="test_type,circuit_A,circuit_B,qubits,ops_A,ops_B"
for s in "${STRATEGIES[@]}"; do
    CSV_HEADER="${CSV_HEADER},${s}_result,${s}_time_s,${s}_nodes"
done
echo "$CSV_HEADER" > "$RESULTS_FILE"

# Run our equiv_check, return result,time,sifted_nodes
run_ours() {
    local file1="$1"
    local file2="$2"
    local strat="$3"
    local output
    output=$(timeout 60 "$EQUIV_CHECK" "$file1" "$file2" "$strat" 2>/dev/null)

    if [ $? -eq 124 ]; then
        echo "TIMEOUT,-1,-1"
        return
    fi

    local result="UNKNOWN"
    if echo "$output" | grep -q "NOT EQUIVALENT"; then
        result="not_equiv"
    elif echo "$output" | grep -q "EQUIVALENT"; then
        result="equiv"
    fi

    local time_val
    time_val=$(echo "$output" | grep "总时间" | grep -oP '[\d.]+(?=s)')
    [ -z "$time_val" ] && time_val="-1"

    # Extract final node count (after sifting if applied)
    local nodes="-1"
    local sift_line=$(echo "$output" | grep "Sifting 后")
    if [ -n "$sift_line" ]; then
        nodes=$(echo "$sift_line" | grep -oP '\d+(?= nodes)')
    else
        # No sifting - use product node count
        local prod_line=$(echo "$output" | grep "A \* B†")
        [ -n "$prod_line" ] && nodes=$(echo "$prod_line" | grep -oP '\d+(?= nodes)')
    fi

    echo "${result},${time_val},${nodes}"
}

# Get circuit info
get_info() {
    local file="$1"
    local qubits=$(grep -oP '(?<=\.numvars )\d+' "$file" 2>/dev/null || echo "0")
    local ops=$(grep -c "^[a-z]" "$file" 2>/dev/null || echo "0")
    echo "$qubits,$ops"
}

# Print table header
print_header() {
    printf "  %-26s %-26s" "Circuit_A" "Circuit_B"
    for s in "${STRATEGIES[@]}"; do
        printf " | %-22s" "$s"
    done
    printf "\n"
    printf "  %-26s %-26s" "" ""
    for s in "${STRATEGIES[@]}"; do
        printf " | %-7s %-7s %-5s" "result" "time" "nodes"
    done
    printf "\n"
    printf "  "
    printf -- "-%.0s" {1..52}
    for s in "${STRATEGIES[@]}"; do
        printf "-+-"
        printf -- "-%.0s" {1..21}
    done
    printf "\n"
}

# Run a single test across all strategies
run_test() {
    local type="$1"
    local name1="$2"
    local name2="$3"
    local file1="$CIRCUITS_DIR/$name1"
    local file2="$CIRCUITS_DIR/$name2"

    [ ! -f "$file1" ] && { echo "  [SKIP] $name1 not found"; return; }
    [ ! -f "$file2" ] && { echo "  [SKIP] $name2 not found"; return; }

    local info1=$(get_info "$file1")
    local info2=$(get_info "$file2")
    local qubits=$(echo "$info1" | cut -d',' -f1)
    local ops1=$(echo "$info1" | cut -d',' -f2)
    local ops2=$(echo "$info2" | cut -d',' -f2)

    printf "  %-26s %-26s" "$name1" "$name2"

    # CSV row start
    local csv_row="${type},${name1},${name2},${qubits},${ops1},${ops2}"

    # Run each strategy
    for strat in "${STRATEGIES[@]}"; do
        local ours=$(run_ours "$file1" "$file2" "$strat")
        local eq=$(echo "$ours" | cut -d',' -f1)
        local t=$(echo "$ours" | cut -d',' -f2)
        local n=$(echo "$ours" | cut -d',' -f3)
        printf " | %-7s %-7s %-5s" "$eq" "${t}s" "$n"
        csv_row="${csv_row},${eq},${t},${n}"
    done

    printf "\n"
    echo "$csv_row" >> "$RESULTS_FILE"
}

# ===== Test Suite =====
echo "--- 测试1: 自等价性 (A == A) ---"
echo ""
print_header

run_test "self_equiv" "3_17_13.real" "3_17_13.real"
run_test "self_equiv" "rd53_130.real" "rd53_130.real"
run_test "self_equiv" "4mod5-bdd_287.real" "4mod5-bdd_287.real"
run_test "self_equiv" "alu-bdd_288.real" "alu-bdd_288.real"
run_test "self_equiv" "f2_232.real" "f2_232.real"
run_test "self_equiv" "cm82a_208.real" "cm82a_208.real"
run_test "self_equiv" "rd53_251.real" "rd53_251.real"
run_test "self_equiv" "hwb8_114.real" "hwb8_114.real"
run_test "self_equiv" "hwb8_113.real" "hwb8_113.real"
run_test "self_equiv" "urf2_154.real" "urf2_154.real"

echo ""
echo "--- 测试2: 非等价性 (A != B) ---"
echo ""
print_header

run_test "not_equiv" "3_17_13.real" "3_17_14.real"
run_test "not_equiv" "rd53_130.real" "alu-bdd_288.real"
run_test "not_equiv" "f2_232.real" "cm82a_208.real"
run_test "not_equiv" "hwb8_114.real" "hwb8_113.real"
run_test "not_equiv" "rd53_138.real" "rd53_139.real"

echo ""
echo "=============================================="
echo " 完成! 结果保存至: $RESULTS_FILE"
echo "=============================================="
echo ""
echo "CSV 内容:"
cat "$RESULTS_FILE"
