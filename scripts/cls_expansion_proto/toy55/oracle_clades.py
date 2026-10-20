#!/usr/bin/env python3
"""Best possible per-particle purity: per-particle scores against the TRUE conformer maps at the refined poses (arm_llr.py
npz: SSE = -log likelihood up to a constant), assigned to the best conformer (MAP) or the best clade (log-sum-exp of the
clade's conformers); accuracy at conformer level and at reference-clade cuts, single images and pooled same-conformer
groups. usage: oracle_clades.py ARMLLR_NPZ TRUTH_NPZ [CUTS]"""
import sys, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
z = np.load(sys.argv[1]); t = np.load(sys.argv[2]); cuts = [float(c) for c in (sys.argv[3] if len(sys.argv) > 3 else '0.45,0.40,0.35,0.30').split(',')]
S = z['SSE'].astype(np.float64); kt = z['ktrue']; confs = list(z['confs']); n = len(kt); rng = np.random.default_rng(0)
Z = linkage(squareform(t['D'], checks=False), 'average'); tids = list(t['ids'])
print('%d particles x %d conformers; conformer-level MAP accuracy %.3f (chance %.3f)' % (n, len(confs), np.mean(np.argmin(S, 1) == kt), 1. / len(confs)))
for cut in cuts:
    cl = fcluster(Z, t=cut, criterion='distance'); cof = np.array([cl[tids.index(k)] for k in confs]); nc = cl.max()
    tc = cof[kt]; pc_map = cof[np.argmin(S, 1)]
    L = -S; lse = np.array([np.logaddexp.reduce(L[:, cof == c], 1) - np.log((cof == c).sum()) if (cof == c).any() else np.full(n, -np.inf) for c in range(1, nc + 1)]).T; pc_lse = np.argmax(lse, 1) + 1   # uniform within the clade: no bonus for clades with more conformers
    share = np.bincount(tc, minlength=nc + 1)[1:] / n
    pooled = []
    for m in (10, 100):
        hit = []
        for c in np.unique(tc):
            rows = np.where(tc == c)[0]
            for _ in range(20): g = rng.choice(rows, m, replace=True); hit.append(int(np.argmax(lse[g].sum(0)) + 1 == c))
        pooled.append(np.mean(hit))
    print('clades at %.2f: %2d clades (largest share %.2f) | per-image accuracy MAP %.3f, clade log-sum-exp %.3f | pooled 10 images %.2f, 100 images %.2f' % (cut, nc, share.max(), np.mean(pc_map == tc), np.mean(pc_lse == tc), pooled[0], pooled[1]))
