#!/usr/bin/env python3
"""Truth readout of a states_clusters.txt on a toy: per state the majority conformer, its purity and size; coverage =
conformers that are the majority of some state; junk = conformers absent from gt (0). usage: score_states.py STATES GT [--tag T]"""
import sys, numpy as np
st = np.loadtxt(sys.argv[1], dtype=int); rows = [l.split() for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
n = max(len(st), max(int(w[0]) for w in rows)); gt = np.zeros(n, int); st = np.r_[st, np.zeros(n - len(st), int)]
for w in rows: gt[int(w[0]) - 1] = int(w[1])
tag = sys.argv[4] if len(sys.argv) > 4 else sys.argv[1]; cov = set(); rows = []
for s in sorted(set(st) - {0}):
    g = gt[st == s]; vals, cnt = np.unique(g, return_counts=True); k = cnt.argmax(); rows.append((s, len(g), vals[k], cnt[k] / len(g), (g == 0).mean())); cov.add(vals[k]) if vals[k] else None
print('%s: %d states; coverage %d conformers %s; mean purity %.2f' % (tag, len(rows), len(cov), sorted(cov), np.mean([r[3] for r in rows]) if rows else 0))
for r in rows: print('  state %d: %d particles, majority conformer %d at purity %.2f, junk %.2f' % r)
