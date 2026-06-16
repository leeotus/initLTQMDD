#!/bin/bash
# Run benchmark_choi_matrix for each (gate, n) combination individually
# to avoid OOM from one combo crashing everything else.
# Results are appended to the same CSV file.
#
# Usage: bash test/run_choi_compare_upto16.sh

BENCH="./build_choi_release/test/benchmark_choi_matrix"
OUT="results/choi_matrix/choi_results_latest.csv"
ERRLOG="results/choi_matrix/choi_errors.log"
TIMEOUT_PER_GATE=300  # 5 minutes per (gate, n)

mkdir -p results/choi_matrix

# If output exists but is missing the header, write it
if [ ! -f "$OUT" ] || [ ! -s "$OUT" ]; then
    echo "Creating new results file with header..."
    $BENCH --nmin 2 --nmax 2 2>/dev/null | head -1 > "$OUT"
fi

GATES=("H_layer" "CNOT_chain" "QFT" "Grover" "Clifford")

for n in $(seq 11 16); do
    for gate in "${GATES[@]}"; do
        if [ "$gate" = "Grover" ] && [ "$n" -lt 2 ]; then continue; fi
        if [ "$gate" = "Clifford" ] && [ "$n" -gt 11 ]; then
            echo "[skip] $gate n=$n (too large, likely OOM)"
            continue
        fi

        echo "[run] $gate n=$n ..."
        result=$(timeout $TIMEOUT_PER_GATE $BENCH --nmin "$n" --nmax "$n" 2>>"$ERRLOG" | grep "^${gate},")
        if [ $? -eq 0 ] && [ -n "$result" ]; then
            echo "$result" >> "$OUT"
            echo "  ✓ done (line appended)"
        else
            echo "  ✗ timeout or error (see $ERRLOG)"
        fi
    done
done

echo ""
echo "=== Done. Total lines: $(wc -l < "$OUT") ==="
echo "Last 5 lines:"
tail -5 "$OUT"
