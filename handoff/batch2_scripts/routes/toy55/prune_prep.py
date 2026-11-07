#!/usr/bin/env python3
"""Pruning test input for one f692e board (region or zoom): its per-meta-row state labels -> per-half candidate labels
(half rows), sub-states with fewer than FLOOR particles in either half left out (too little signal to test), one
numbering for both halves. usage: prune_prep.py LABELS META_ROWS_A META_ROWS_B FLOOR OUTDIR"""
import sys, os, numpy as np
labf, mra, mrb, floor, out = sys.argv[1:6]; floor = int(floor); os.makedirs(out, exist_ok=True)
L = np.loadtxt(labf).astype(int); half = {h: L[np.loadtxt(f).astype(int) - 1] for h, f in (('A', mra), ('B', mrb))}
ids = sorted(set(L[L > 0])); cnt = {k: (int((half['A'] == k).sum()), int((half['B'] == k).sum())) for k in ids}
keep = [k for k in ids if min(cnt[k]) >= floor]; small = [k for k in ids if k not in keep]
for h in 'AB':
    c = np.zeros(len(half[h]), int)
    for i, k in enumerate(keep): c[half[h] == k] = i + 1
    np.savetxt('%s/cand_%s.txt' % (out, h), c, fmt='%d')
with open(out + '/labels.txt', 'w') as f:
    f.write('# cand original nA nB\n'); [f.write('%d %d %d %d\n' % (i + 1, k, *cnt[k])) for i, k in enumerate(keep)]
with open(out + '/small.txt', 'w') as f:
    f.write('# original nA nB (below the floor of %d per half)\n' % floor); [f.write('%d %d %d\n' % (k, *cnt[k])) for k in small]
print('%d sub-states: %d tested, %d below %d per half' % (len(ids), len(keep), len(small), floor))
