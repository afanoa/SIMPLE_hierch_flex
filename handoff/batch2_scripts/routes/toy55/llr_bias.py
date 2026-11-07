#!/usr/bin/env python3
"""Is a chance-level arm_llr result a per-map bias? Scores column-centred (per-conformer mean over all images removed),
then accuracy, truth rank and pooled accuracy again. usage: llr_bias.py NPZ"""
import sys, numpy as np
z = np.load(sys.argv[1]); kt = z['ktrue']; rng = np.random.default_rng(0)
for tag in ('SSE', 'SSEm'):
    S0 = z[tag]; n, K = S0.shape
    for name, S in (('raw', S0), ('column-centred', S0 - S0.mean(0))):
        top = np.bincount(np.argmin(S, 1), minlength=K)
        print('%s %-15s accuracy %.3f (chance %.3f), truth rank median %d, most-chosen map takes %.2f' % (tag, name, np.mean(np.argmin(S, 1) == kt), 1. / K, np.median([(S[i] < S[i, kt[i]]).sum() + 1 for i in range(n)]), top.max() / n))
        for m in (10, 100, 1000):
            hit = []
            for t in range(K):
                rows = np.where(kt == t)[0]
                if len(rows) < 8: continue
                for _ in range(20): g = rng.choice(rows, m, replace=True); hit.append(int(np.argmin(S[g].sum(0)) == t))
            print('    pooled %4d: accuracy %.2f' % (m, np.mean(hit)))
