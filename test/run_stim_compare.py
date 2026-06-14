#!/usr/bin/env python3
"""Quick Stim comparison: Steane [[7,1,3]] QEC circuits, noise-free."""
import time, csv, os, stim
from datetime import datetime

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'qec_comparison_results')
os.makedirs(RESULTS_DIR, exist_ok=True)

def build_stim_steane(rounds):
    c = stim.Circuit()
    # Encode
    c.append('H', [0]); c.append('H', [1]); c.append('H', [3])
    for ctrl, tgt in [(0,2),(0,4),(1,2),(1,5),(3,4),(3,5),(2,6),(4,6),(5,6)]:
        c.append('CNOT', [ctrl, tgt])
    # Syndrome extraction rounds
    xs = [[0,2,4,6],[1,2,5,6],[3,4,5,6]]
    zs = [[0,2,4,6],[1,2,5,6],[3,4,5,6]]
    for _ in range(rounds):
        for s,qs in enumerate(xs):
            a=7+s; c.append('H',[a])
            for q in qs: c.append('CNOT',[a,q])
            c.append('H',[a])
        for s,qs in enumerate(zs):
            a=10+s
            for q in qs: c.append('CNOT',[q,a])
    for q in range(7): c.append('M',[q])
    return c

results = []
for rounds in [1, 5, 10, 20]:
    c = build_stim_steane(rounds)
    t0 = time.time()
    # Use 100 shots for speed
    sampler = c.compile_sampler()
    data = sampler.sample(100)
    elapsed = time.time() - t0
    errors = sum(1 for row in data if sum(int(row[i]) for i in range(7)) % 2 == 1)
    ler = errors / 100
    results.append({'rounds':rounds, 'stim_time_s':f'{elapsed:.4f}', 'stim_ler':f'{ler:.2f}'})
    print(f'Rounds={rounds}: time={elapsed:.4f}s, LER={ler:.2f}')

ts = datetime.now().strftime('%Y%m%d_%H%M%S')
csv_path = os.path.join(RESULTS_DIR, f'stim_compare_{ts}.csv')
with open(csv_path,'w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=['rounds','stim_time_s','stim_ler'])
    w.writeheader()
    w.writerows(results)
print(f'Saved to {csv_path}')