#!/usr/bin/env python3
"""True-conformer state labels for a tree's root rows: one label per project row (1..K in conformer order, 0 = outside
the root or junk). usage: oracle_states.py GT SUBSET OUT"""
import sys, numpy as np
gt, sub, out = sys.argv[1:4]; conf = np.loadtxt(gt, comments='#')[:, 1].astype(int); rows = np.loadtxt(sub, dtype=int) - 1
ks = sorted(set(conf[rows]) - {0}); st = np.zeros(len(conf), int)
for i, k in enumerate(ks): st[rows[conf[rows] == k]] = i + 1
np.savetxt(out, st, fmt='%d'); print('%d rows, conformers %s -> states 1..%d, sizes %s' % ((st > 0).sum(), ks, len(ks), np.bincount(st)[1:].tolist()))
