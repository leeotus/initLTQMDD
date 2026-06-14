#!/bin/bash
# QEC Benchmark Script
# Runs all 5 experiments using the QMDD-based QEC simulation tool

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BIN="${PROJECT_DIR}/build/qec_benchmark"
OUTPUT_DIR="${PROJECT_DIR}/qec_results"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=============================================="
echo "  QEC Benchmark - QMDD-based QEC Simulation"
echo "  Steane [[7,1,3]] + QMDD + ddsim noise engine"
echo "=============================================="
echo ""

# Check binary exists
if [ ! -f "$BIN" ]; then
    echo -e "${RED}[ERROR] qec_benchmark not found at $BIN${NC}"
    echo ""
    echo "Build it with:"
    echo "  cd $PROJECT_DIR"
    echo "  cmake -S . -B build_qec -DBUILD_QEC=ON -DGIT_SUBMODULE=OFF -DCMAKE_BUILD_TYPE=Release"
    echo "  cmake --build build_qec --target qec_benchmark -j\$(nproc)"
    echo "  cp build_qec/qec/qec_benchmark build/qec_benchmark"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${OUTPUT_DIR}/qec_benchmark_${TIMESTAMP}.log"

echo "Date:       $(date)"
echo "Binary:     $BIN"
echo "Output:     $LOG_FILE"
echo ""

# Run all 5 experiments
echo -e "${GREEN}[1/5] Running Experiment 1: Encoding Verification...${NC}"
echo -e "${GREEN}[2/5] Running Experiment 2: DD Compression...${NC}"
echo -e "${GREEN}[3/5] Running Experiment 3: Noise Threshold Sweep...${NC}"
echo -e "${GREEN}[4/5] Running Experiment 4: DD Scalability vs QEC Rounds...${NC}"
echo -e "${GREEN}[5/5] Running Experiment 5: IG Symmetry Analysis...${NC}"
echo ""

$BIN 2>&1 | tee "$LOG_FILE"

echo ""
echo "=============================================="
echo "  Benchmark completed!"
echo "  Full log: $LOG_FILE"
echo "=============================================="

# Extract summary statistics
echo ""
echo -e "${YELLOW}Summary Statistics:${NC}"
echo ""

# Experiment 1: PASS/FAIL
if grep -q "PASSED" "$LOG_FILE"; then
    echo -e "  Experiment 1 (Encoding):  ${GREEN}PASSED${NC}"
else
    echo -e "  Experiment 1 (Encoding):  ${RED}FAILED${NC} (check circuit implementation)"
fi

# Experiment 2: Compression ratio
COMP_RATIO=$(grep "Compression ratio:" "$LOG_FILE" | awk '{print $NF}' | sed 's/x//')
echo "  Experiment 2 (DD Comp):  Compression ratio = ${COMP_RATIO}"

# Experiment 3: Noise threshold data
echo "  Experiment 3 (Noise):    Data points collected:"
grep -E "p_physical|^[0-9]" "$LOG_FILE" | sed 's/^/    /'

# Experiment 4: Rounds data
echo "  Experiment 4 (Rounds):   Data points collected:"
grep -E "Rounds|^[0-9]" "$LOG_FILE" | sed 's/^/    /'

# Experiment 5: Symmetry groups
SYM_COUNT=$(grep "Symmetry groups detected:" "$LOG_FILE" | awk '{print $NF}')
echo "  Experiment 5 (IG Sym):   Symmetry groups detected = ${SYM_COUNT}"

echo ""
echo "=============================================="
echo "  Done"
echo "=============================================="