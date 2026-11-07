#!/usr/bin/env python3
"""usage: full692_rank_check.py TABLE REF_LABELS STATE: every site pole's posterior of a factorial table ranked against a reference particle
set (REF_LABELS == STATE, e.g. the f692st closed leaf): AUC over the pole's level and the reference share among its top 2.5k / 10k."""
import sys, numpy as np
TAB, RF, ST = sys.argv[1], sys.argv[2], int(sys.argv[3]); ref = np.loadtxt(RF, dtype=int) == ST; NM = len(ref)
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
R = []
for w in [l.split() for l in open(TAB) if l.strip() and not l.startswith('#')]:
    name, lev, labf, postf, rowf = w[:5]; Q = np.loadtxt(postf, ndmin=2); Q = Q.T if Q.shape[0] != NM else Q
    r = np.loadtxt(rowf, dtype=int).ravel(); rows = (r > 0) if (len(r) == NM and r.max() <= 1) else np.isin(np.arange(NM), r - 1)
    for j in range(Q.shape[1]):
        P = Q[:, j]; ok = rows & (P >= 0)
        if ref[ok].sum() < 50: continue
        a = auc(P[ok], ref[ok]); s = P if a >= 0.5 else -P; o = np.argsort(-np.where(ok, s, -np.inf))
        R.append((max(a, 1 - a), name, j + 1, lev, ok.sum(), ref[ok].sum(), ref[o[:2500]].mean(), ref[o[:10000]].mean(), ref[ok].mean()))
print('%d site poles ranked against the reference set (%d particles)' % (len(R), ref.sum()))
for r in sorted(R, reverse=True)[:12]: print('  %-10s pole %d level %-7s AUC %.3f | top-2.5k %.3f top-10k %.3f (base %.3f) | %d of the reference in the level' % (r[1], r[2], r[3], r[0], r[6], r[7], r[8], r[5]))
