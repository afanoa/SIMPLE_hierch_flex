#!/usr/bin/env python3
"""Level-2 ranking without truth: the cells of a within-set board scored by the correlation of their (half-averaged)
deviation with a TARGET deviation (the parent cluster's, from cluster_targets.npz), every particle by its
responsibility-weighted cell score; the top N get state 1 in OUT/states_top.txt; the per-particle scores saved.
usage: rank_by_target.py BOARD TARGETS_NPZ --target k --n N --nrows NROWS --out OUT"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('targets', nargs='?', default=''); ap.add_argument('--target', default='1', help='parent state number, or pc1: the ends of the within-board deviations first principal axis (truth-free, no parent needed)'); ap.add_argument('--n', type=int, default=10000); ap.add_argument('--nrows', type=int, required=True); ap.add_argument('--out', required=True); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); H = np.load(a.board + '/klem_halves.npz'); D = 0.5 * (H['D0'] + H['D1'])
if a.target == 'pc1':   # the board's own main deviation axis: first right-singular vector of the mean-centred cell deviations (cross-half symmetrised)
    Dc = D - D.mean(0, keepdims=True); U, S_, Vt = np.linalg.svd(Dc, full_matrices=False); t = Vt[0]
    cs = np.real(Dc @ np.conj(t)); cs -= cs.mean(); print('own main axis: singular values (first 4, relative) %s' % np.array2string(S_[:4] / S_[0], precision=2))
else:
    Tz = np.load(a.targets); k = list(Tz['states']).index(int(a.target)); t = Tz['T'][k]
    assert int(Tz['n3']) == int(H['n3']) and np.array_equal(Tz['band'], H['band']), 'band mismatch between the parent and the within-set board'
    cs = np.array([np.real(np.sum(np.conj(d) * t)) / max(np.sqrt(np.sum(np.abs(d) ** 2) * np.sum(np.abs(t) ** 2)), 1e-30) for d in D]); cs -= cs.mean()
R = np.load(a.board + '/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(D); s = r[:, :ngrid] @ cs; order = np.argsort(-s)
st = np.zeros(a.nrows, int); st[pinds[order[:a.n]] - 1] = 1
if a.target == 'pc1': st[pinds[order[-a.n:]] - 1] = 2   # both ends of the axis
np.savetxt(a.out + '/states_top.txt', st, fmt='%d'); np.savetxt(a.out + '/score_per_particle.txt', np.c_[pinds, s], fmt='%d %.4f')
print('cells scored against %s: cell score range %+.2f .. %+.2f; particle score quantiles 10/50/90 = %+.3f %+.3f %+.3f; top %d >= %+.3f' % ('the own main axis (pc1, both ends)' if a.target == 'pc1' else 'the parent state ' + a.target, cs.min(), cs.max(), *np.percentile(s, [10, 50, 90]), a.n, s[order[a.n - 1]]))
