#!/usr/bin/env python3
"""
QEC Comparison Script: Ours (QMDD+QEC) vs Stim vs Qiskit Aer
Compares simulation time, LER, and DD-specific metrics on Steane [[7,1,3]] code.
"""

import subprocess, time, csv, os, sys, json, math
from datetime import datetime
import numpy as np

# ================== Configuration ==================
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUR_BIN = os.path.join(PROJECT_DIR, "build", "qec_benchmark")
RESULTS_DIR = os.path.join(PROJECT_DIR, "qec_comparison_results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Check tool availability
try:
    import stim
    HAS_STIM = True
    print(f"stim: {stim.__version__}")
except ImportError:
    HAS_STIM = False
    print("stim: NOT INSTALLED")

try:
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator
    HAS_AER = True
    print("qiskit-aer: OK")
except ImportError:
    HAS_AER = False
    print("qiskit-aer: NOT INSTALLED")

print(f"Ours (QMDD): {OUR_BIN} — {'OK' if os.path.isfile(OUR_BIN) else 'NOT FOUND'}")
print()
if not os.path.isfile(OUR_BIN):
    print("ERROR: qec_benchmark not found. Build it first.")
    sys.exit(1)

# ================== Steane [[7,1,3]] Circuit Definition ==================
# Steane code: 7 data qubits + 6 ancilla qubits = 13 total
# X-stabilizers (from classical Hamming [7,4,3]):
#   S0: X0 X2 X4 X6
#   S1: X1 X2 X5 X6
#   S2: X3 X4 X5 X6
# Z-stabilizers:
#   S3: Z0 Z2 Z4 Z6
#   S4: Z1 Z2 Z5 Z6
#   S5: Z3 Z4 Z5 Z6

def build_steane_qec_circuit_stim(num_rounds, noise_p=0.0):
    """Build Steane [[7,1,3]] QEC circuit using Stim."""
    c = stim.Circuit()
    n_data = 7
    n_anc = 6
    n_total = n_data + n_anc

    # Step 1: Encode logical |0⟩_L
    # H on qubits 0,1,3; then CNOTs for Hamming encoder
    c.append("H", [0])
    c.append("H", [1])
    c.append("H", [3])
    # CNOT gates
    cnots = [(0,2),(0,4),(1,2),(1,5),(3,4),(3,5),(2,6),(4,6),(5,6)]
    for ctrl, tgt in cnots:
        c.append("CNOT", [ctrl, tgt])

    # Step 2: QEC rounds (syndrome extraction)
    x_stabs = [
        [0,2,4,6],
        [1,2,5,6],
        [3,4,5,6],
    ]
    z_stabs = [
        [0,2,4,6],
        [1,2,5,6],
        [3,4,5,6],
    ]

    for rnd in range(num_rounds):
        # X-type stabilizer measurements
        for s_idx, qubits in enumerate(x_stabs):
            anc = n_data + s_idx  # ancilla 0-2 for X
            c.append("H", [anc])
            for q in qubits:
                c.append("CNOT", [anc, q])
            c.append("H", [anc])

            # Noise injection
            if noise_p > 0:
                c.append("DEPOLARIZE1", [anc], noise_p)

        # Z-type stabilizer measurements
        for s_idx, qubits in enumerate(z_stabs):
            anc = n_data + 3 + s_idx  # ancilla 3-5 for Z
            for q in qubits:
                c.append("CNOT", [q, anc])

            # Noise injection
            if noise_p > 0:
                c.append("DEPOLARIZE1", [anc], noise_p)

        # Depolarizing noise on each data qubit per round
        if noise_p > 0:
            for q in range(n_data):
                c.append("DEPOLARIZE1", [q], noise_p)

    # Step 3: Measure all data qubits in Z basis
    for q in range(n_data):
        c.append("M", [q])

    return c

# ================== Stim Simulation ==================

def run_stim_logical_errors(circuit, shots=1000):
    """Run Stim sampler and compute logical error rate for Steane code.
    Logical Z parity = XOR of all 7 data qubit measurements.
    """
    sampler = circuit.compile_sampler()
    results = sampler.sample(shots)

    logical_errors = 0
    for row in results:
        # row = measurements of [q0, q1, ..., q6, anc0, anc1, ...]  (ancilla quantities depend on M ops)
        # We only have 7 data qubits measured (the M operations)
        parity = 0
        for i in range(7):
            parity ^= int(row[i])
        if parity == 1:  # logical |1> instead of logical |0>
            logical_errors += 1

    ler = logical_errors / shots
    return ler

# ================== Qiskit Aer Simulation ==================

def build_steane_qec_circuit_aer(num_rounds, noise_p=0.0):
    """Build Steane [[7,1,3]] QEC circuit using Qiskit."""
    from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
    n_data = 7
    n_anc = 6
    n_total = n_data + n_anc

    qr = QuantumRegister(n_total, 'q')
    cr = ClassicalRegister(n_data, 'c')
    qc = QuantumCircuit(qr, cr)

    # Encode
    qc.h(0); qc.h(1); qc.h(3)
    cnots = [(0,2),(0,4),(1,2),(1,5),(3,4),(3,5),(2,6),(4,6),(5,6)]
    for ctrl, tgt in cnots:
        qc.cx(ctrl, tgt)

    x_stabs = [[0,2,4,6],[1,2,5,6],[3,4,5,6]]
    z_stabs = [[0,2,4,6],[1,2,5,6],[3,4,5,6]]

    noise_model = None
    if noise_p > 0:
        from qiskit_aer.noise import NoiseModel, depolarizing_error
        noise_model = NoiseModel()
        error_1q = depolarizing_error(noise_p, 1)
        error_2q = depolarizing_error(noise_p, 2)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'x', 'z', 'measure'])
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx'])
        noise_model.add_all_qubit_readout_error(
            [1-noise_p, noise_p])

    for rnd in range(num_rounds):
        for s_idx, qubits in enumerate(x_stabs):
            anc = n_data + s_idx
            qc.h(anc)
            for q in qubits:
                qc.cx(anc, q)
            qc.h(anc)

        for s_idx, qubits in enumerate(z_stabs):
            anc = n_data + 3 + s_idx
            for q in qubits:
                qc.cx(q, anc)

    # Measure all data qubits
    for i in range(n_data):
        qc.measure(qr[i], cr[i])

    if noise_p > 0:
        return qc, noise_model
    return qc, None

# ================== Ours (QMDD) Simulation ==================

def run_ours(num_rounds):
    """Run our QMDD+QEC benchmark and parse output."""
    # We modify qec_benchmark to support command-line args for rounds
    # For now, use a temp approach: create a wrapper
    t0 = time.time()
    result = subprocess.run(
        [OUR_BIN],
        capture_output=True, text=True, timeout=300,
        cwd=PROJECT_DIR
    )
    elapsed = time.time() - t0

    # Parse output for key metrics
    output = result.stdout

    # Experiment 1: PASS/FAIL
    passed = "PASSED" in output

    # Experiment 2: Compression ratio
    comp_ratio = 1.0
    dd_size_no_sift = 0
    dd_size_sift = 0
    for line in output.split('\n'):
        if "DD Size (no sifting):" in line:
            parts = line.split()
            dd_size_no_sift = int(parts[4])
        if "DD Size (with sifting):" in line:
            parts = line.split()
            dd_size_sift = int(parts[4])
        if "Compression ratio:" in line:
            parts = line.split()
            comp_ratio = float(parts[2].replace('x', ''))

    # Experiment 3: LER data (last line = highest noise level)
    ler_data = []
    in_table = False
    for line in output.split('\n'):
        if "p_physical" in line:
            in_table = True
            continue
        if in_table and line.strip() and line[0].isdigit():
            parts = line.split()
            if len(parts) >= 5:
                try:
                    p = float(parts[0])
                    ler = float(parts[1])
                    ler_data.append((p, ler))
                except ValueError:
                    pass

    return {
        "time": elapsed,
        "passed": passed,
        "comp_ratio": comp_ratio,
        "dd_size_no_sift": dd_size_no_sift,
        "dd_size_sift": dd_size_sift,
        "ler_data": ler_data,
    }

# ================== Main Comparison Loop ==================

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = os.path.join(RESULTS_DIR, f"comparison_{timestamp}.csv")

    # Experiment configurations
    configs = [
        {"rounds": 1, "noise": 0.0, "shots": 2000, "label": "1 round, noise-free"},
        {"rounds": 5, "noise": 0.0, "shots": 2000, "label": "5 rounds, noise-free"},
        {"rounds": 10, "noise": 0.0, "shots": 2000, "label": "10 rounds, noise-free"},
        {"rounds": 1, "noise": 0.001, "shots": 1000, "label": "1 round, p=0.001"},
        {"rounds": 5, "noise": 0.001, "shots": 1000, "label": "5 rounds, p=0.001"},
        {"rounds": 5, "noise": 0.01, "shots": 1000, "label": "5 rounds, p=0.01"},
    ]

    results = []

    print("=" * 80)
    print("  QEC Comparison: Ours (QMDD+QEC) vs Stim vs Qiskit Aer")
    print(f"  Circuit: Steane [[7,1,3]] (13 physical qubits)")
    print(f"  Output:  {csv_path}")
    print("=" * 80)
    print()

    for cfg in configs:
        print(f"--- {cfg['label']} ---")

        row = {
            "rounds": cfg["rounds"],
            "noise_p": cfg["noise"],
            "shots": cfg["shots"],
        }

        # ===== Ours (QMDD) =====
        print("  Ours (QMDD)...", end=" ", flush=True)
        t0 = time.time()
        our = run_ours(cfg["rounds"])
        row["ours_time"] = our["time"]
        row["ours_dd_size"] = our["dd_size_sift"] if our["dd_size_sift"] > 0 else our["dd_size_no_sift"]
        row["ours_comp_ratio"] = our["comp_ratio"]

        # Get LER at closest noise level
        ours_ler = None
        for p, ler in our["ler_data"]:
            if abs(p - cfg["noise"]) < 1e-6:
                ours_ler = ler
                break
        row["ours_ler"] = ours_ler if ours_ler is not None else "N/A"
        print(f"{our['time']:.2f}s, LER={ours_ler}")

        # ===== Stim =====
        if HAS_STIM:
            print("  Stim...    ", end=" ", flush=True)
            t0 = time.time()
            c = build_steane_qec_circuit_stim(cfg["rounds"], cfg["noise"])
            ler = run_stim_logical_errors(c, cfg["shots"])
            stim_time = time.time() - t0
            row["stim_time"] = stim_time
            row["stim_ler"] = ler
            print(f"{stim_time:.3f}s, LER={ler:.4f}")
        else:
            row["stim_time"] = "N/A"
            row["stim_ler"] = "N/A"
            print("  Stim...    SKIP (not installed)")

        # ===== Qiskit Aer =====
        if HAS_AER and cfg["rounds"] <= 5:
            print("  Aer...     ", end=" ", flush=True)
            t0 = time.time()
            qc, noise_model = build_steane_qec_circuit_aer(cfg["rounds"], cfg["noise"])
            backend = AerSimulator(method='automatic')
            if noise_model:
                job = backend.run(qc, noise_model=noise_model, shots=cfg["shots"])
            else:
                job = backend.run(qc, shots=cfg["shots"])
            counts = job.result().get_counts()
            aer_time = time.time() - t0

            # Compute LER from counts
            logical_errors = 0
            for bitstr, cnt in counts.items():
                # bitstr is reversed in qiskit (q0 is rightmost)
                parity = 0
                bitstr_rev = bitstr[::-1]  # reverse to get q0 first
                for i in range(7):
                    if i < len(bitstr_rev) and bitstr_rev[i] == '1':
                        parity ^= 1
                if parity == 1:
                    logical_errors += cnt
            aer_ler = logical_errors / cfg["shots"]
            row["aer_time"] = aer_time
            row["aer_ler"] = aer_ler
            print(f"{aer_time:.2f}s, LER={aer_ler:.4f}")
        else:
            row["aer_time"] = "N/A" if not HAS_AER else "SKIP (too many rounds)"
            row["aer_ler"] = "N/A" if not HAS_AER else "SKIP"
            print("  Aer...     SKIP")

        results.append(row)
        print()

    # Write CSV
    fieldnames = [
        "rounds", "noise_p", "shots",
        "ours_time", "ours_ler", "ours_dd_size", "ours_comp_ratio",
        "stim_time", "stim_ler",
        "aer_time", "aer_ler",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print("=" * 80)
    print(f"  Results saved to {csv_path}")
    print("=" * 80)
    print()

    # Print summary table
    print(f"{'Rounds':<8} {'Noise':<10} {'Ours(s)':<10} {'Ours LER':<12} {'Stim(s)':<10} {'Stim LER':<12} {'Aer(s)':<10} {'Aer LER':<12}")
    print("-" * 90)
    for r in results:
        print(f"{r['rounds']:<8} {r['noise_p']:<10.0e} {str(r['ours_time']):<10} {str(r['ours_ler']):<12} {str(r['stim_time']):<10} {str(r['stim_ler']):<12} {str(r['aer_time']):<10} {str(r['aer_ler']):<12}")

if __name__ == "__main__":
    main()