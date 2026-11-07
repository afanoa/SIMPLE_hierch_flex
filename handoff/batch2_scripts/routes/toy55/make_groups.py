#!/usr/bin/env python3
"""usage: make_groups.py NROWS OUT SPEC_A SPEC_B
One value per project row: 1 = group A, 2 = group B, 0 = neither. A SPEC is one of
  bins:FILE:K[:ROWS]  rows whose value in FILE is K; FILE holds one value per project row, or per node row with ROWS (the node's
                      rows.txt, 1-based project rows) mapping it
  rows:FILE           the 1-based project rows listed in FILE
  eq:FILE:K           as bins without ROWS (e.g. hetero_tree_labels.txt == leaf index K)"""
import sys, numpy as np
NM, OUT = int(sys.argv[1]), sys.argv[2]; g = np.zeros(NM, int)
def rows_of(spec):
    w = spec.split(':')
    if w[0] == 'rows': return np.loadtxt(w[1], dtype=int).ravel() - 1
    v = np.loadtxt(w[1], dtype=int).ravel(); k = int(w[2])
    if len(w) > 3: r = np.loadtxt(w[3], dtype=int).ravel() - 1; assert len(r) == len(v), '%s: %d values, %d rows' % (w[1], len(v), len(r)); return r[v == k]
    assert len(v) == NM, '%s: %d values, %d project rows' % (w[1], len(v), NM); return np.nonzero(v == k)[0]
a, b = rows_of(sys.argv[3]), rows_of(sys.argv[4]); g[a] = 1; g[b] = 2; g[np.intersect1d(a, b)] = 0
np.savetxt(OUT, g, fmt='%d'); print('%s: A %d, B %d rows' % (OUT, (g == 1).sum(), (g == 2).sum()))
