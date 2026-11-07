#!/usr/bin/env python3
"""usage: rank_agree.py ROWS NAME=POST_MASTER [NAME=POST_MASTER ...]: pairwise Spearman and top-10 % / top-20 % Jaccard over the rows of ROWS (0/1 per master row)"""
import sys, numpy as np
m = np.loadtxt(sys.argv[1], dtype=int) > 0; E = [a.split('=', 1) for a in sys.argv[2:]]; P = {n: np.loadtxt(f) for n, f in E}; ok = m & np.all([p >= 0 for p in P.values()], 0); n = ok.sum()
R = {k: np.argsort(np.argsort(-p[ok])) for k, p in P.items()}; nm = list(P)
print('%d rows (in every ranking)' % n)
for i in range(len(nm)):
    for j in range(i + 1, len(nm)):
        a, b = nm[i], nm[j]; jac = []
        for f in (0.1, 0.2):
            k = int(f * n); A, B = set(np.flatnonzero(R[a] < k)), set(np.flatnonzero(R[b] < k)); jac.append(len(A & B) / len(A | B))
        print('%-14s vs %-14s Spearman %+.3f | Jaccard top-10%% %.2f top-20%% %.2f' % (a, b, np.corrcoef(R[a], R[b])[0, 1], *jac))
