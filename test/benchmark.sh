#!/bin/bash
#
# QMDD Sifting Benchmark Script
#
# Usage:
#   ./benchmark.sh -i <circuit_file_or_dir> -s <strategy1,strategy2,...> -o <output.csv> [-t timeout]
#
# Examples:
#   ./benchmark.sh -i ~/workshop/circuits/ham15_107.real -s sifting,lb,iglb -o results.csv
#   ./benchmark.sh -i ~/workshop/circuits/ -s sifting,lb,iglb,group,iggroup -o bench.csv -t 120
#
# Available strategies:
#   none, sifting, lb, tightlb, ig, iglb, group, iggroup,
#   upperls, lowerls, mixls, lbupperls, lblowerls, lbmixls,
#   igupperls, iglowerls, igmixls

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LTQMDD="$PROJECT_ROOT/build_test/ltqmdd"

TIMEOUT=300
INPUT=""
STRATEGIES=""
OUTPUT=""

usage() {
    echo "Usage: $0 -i <circuit_file_or_dir> -s <strategies> -o <output.csv> [-t timeout_seconds]"
    echo ""
    echo "Options:"
    echo "  -i    Input circuit file (.real/.qasm) or directory containing circuit files"
    echo "  -s    Comma-separated list of sifting strategies"
    echo "  -o    Output CSV file path"
    echo "  -t    Timeout per run in seconds (default: 300)"
    echo ""
    echo "Strategies: none, sifting, lb, tightlb, ig, iglb, group, iggroup,"
    echo "            upperls, lowerls, mixls, lbupperls, lblowerls, lbmixls,"
    echo "            igupperls, iglowerls, igmixls"
    exit 1
}

while getopts "i:s:o:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        s) STRATEGIES="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) TIMEOUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$STRATEGIES" || -z "$OUTPUT" ]]; then
    echo "Error: -i, -s, -o are all required."
    usage
fi

if [[ ! -f "$LTQMDD" ]]; then
    echo "Error: ltqmdd binary not found at $LTQMDD"
    echo "Please build the project first: cd build_test && make ltqmdd"
    exit 1
fi

# Collect circuit files
CIRCUITS=()
if [[ -d "$INPUT" ]]; then
    while IFS= read -r f; do
        CIRCUITS+=("$f")
    done < <(find "$INPUT" -maxdepth 1 -type f \( -name "*.real" -o -name "*.qasm" \) | sort)
elif [[ -f "$INPUT" ]]; then
    CIRCUITS=("$INPUT")
else
    echo "Error: $INPUT is not a valid file or directory."
    exit 1
fi

if [[ ${#CIRCUITS[@]} -eq 0 ]]; then
    echo "Error: No circuit files found in $INPUT"
    exit 1
fi

# Parse strategies
IFS=',' read -ra STRATS <<< "$STRATEGIES"

# Write CSV header
echo -n "circuit,qubits,gates,initial_size" > "$OUTPUT"
for strat in "${STRATS[@]}"; do
    echo -n ",${strat}_size,${strat}_time" >> "$OUTPUT"
done
echo "" >> "$OUTPUT"

echo "=========================================="
echo " QMDD Sifting Benchmark"
echo "=========================================="
echo " Circuits: ${#CIRCUITS[@]} files"
echo " Strategies: ${STRATS[*]}"
echo " Timeout: ${TIMEOUT}s"
echo " Output: $OUTPUT"
echo "=========================================="
echo ""

TOTAL=${#CIRCUITS[@]}
COUNT=0

for circuit in "${CIRCUITS[@]}"; do
    COUNT=$((COUNT + 1))
    filename=$(basename "$circuit")
    name="${filename%.*}"

    printf "[%d/%d] %s ... " "$COUNT" "$TOTAL" "$name"

    # Get initial info with 'none' strategy
    none_out=$(timeout "$TIMEOUT" "$LTQMDD" "$circuit" none 2>&1)
    if [[ $? -ne 0 ]]; then
        echo "SKIP (build failed)"
        continue
    fi

    initial_size=$(echo "$none_out" | grep -oP '初始DD大小: \K[0-9]+')
    qubits=$(echo "$none_out" | grep -oP '(\d+)(?= qubit)' || echo "?")

    # Try to extract qubit/gate count from the .real file header
    if [[ -f "$circuit" ]]; then
        qubits_line=$(grep -i "^\.numvars" "$circuit" 2>/dev/null | grep -oP '\d+')
        if [[ -n "$qubits_line" ]]; then
            qubits="$qubits_line"
        fi
        gates_line=$(grep -c "^t[0-9]\|^f[0-9]\|^T[0-9]\|^F[0-9]" "$circuit" 2>/dev/null)
        if [[ -z "$gates_line" || "$gates_line" == "0" ]]; then
            gates_line=$(grep -c "^[a-z]" "$circuit" 2>/dev/null)
        fi
        gates="${gates_line:-?}"
    else
        qubits="?"
        gates="?"
    fi

    if [[ -z "$initial_size" ]]; then
        initial_size="?"
    fi

    # Write row start
    echo -n "$name,$qubits,$gates,$initial_size" >> "$OUTPUT"

    # Run each strategy
    for strat in "${STRATS[@]}"; do
        strat_out=$(timeout "$TIMEOUT" "$LTQMDD" "$circuit" "$strat" 2>&1)
        exit_code=$?

        if [[ $exit_code -eq 124 ]]; then
            # Timeout
            size="TIMEOUT"
            time_val=">${TIMEOUT}"
        elif [[ $exit_code -ne 0 ]]; then
            # Crash/error
            size="FAIL"
            time_val="FAIL"
        else
            size=$(echo "$strat_out" | grep -oP '大小: \K[0-9]+' | tail -1)
            time_val=$(echo "$strat_out" | grep -oP '时间: \K[0-9.]+' | tail -1)
            if [[ -z "$size" ]]; then
                size="?"
            fi
            if [[ -z "$time_val" ]]; then
                time_val="?"
            fi
        fi

        echo -n ",$size,$time_val" >> "$OUTPUT"
    done

    echo "" >> "$OUTPUT"
    echo "done (init=$initial_size)"
done

echo ""
echo "=========================================="
echo " Benchmark complete. Results in: $OUTPUT"
echo "=========================================="

# Print summary table
echo ""
echo "Summary (DD sizes):"
echo ""
head -1 "$OUTPUT" | tr ',' '\t'
echo "---"
tail -n +2 "$OUTPUT" | while IFS=',' read -r line; do
    echo "$line" | tr ',' '\t'
done | head -20

if [[ $TOTAL -gt 20 ]]; then
    echo "... ($((TOTAL - 20)) more rows in $OUTPUT)"
fi
