#!/usr/bin/env python3
"""lsamp_pair_coords.py: pole-pair coordinates t = s_a - s_b of one run (per pixel set, noise-normalised as lsamp2h.py --ppnorm 1),
optionally over a subset of the node's particles (those owned by the listed states), written as OUT/<tag>/<name>/coords.npz (pole a
low) with an orientation stub. --roles auto names three poles of a K = 4 (or 3) score space c / oA / oB by overlap with reference
posteriors (--ref_post closed,openA,openB); a fourth state is 'rest'. Pair spec TAG:A:B[:S1+S2..] with A / B / S as state names or
the role keywords c, oA, oB, rest.
usage: lsamp_pair_coords.py --scores SC --name RUN --out OUT --pairs TAG:A:B[:S..],.. [--roles auto --ref_post P1,P2,P3]"""
import os, argparse, itertools, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--name', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--pairs', required=True); ap.add_argument('--roles', default=''); ap.add_argument('--ref_post', default='')
a = ap.parse_args()
M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; rows = np.where((OWN >= 0) & (EO >= 0))[0]; K = len(names)
role = {n: n for n in names}
if a.roles == 'auto':
    R = [np.loadtxt(f) for f in a.ref_post.split(',')]; T = np.array([[R[j][rows[OWN[rows] == k]].mean() for j in range(3)] for k in range(K)])
    best = max(itertools.permutations(range(K), 3), key=lambda p: sum(T[p[j], j] for j in range(3))); role.update({'c': names[best[0]], 'oA': names[best[1]], 'oB': names[best[2]]})
    if K > 3: role['rest'] = names[[k for k in range(K) if k not in best][0]]
    print('%s roles: %s' % (a.name, '; '.join('%s %d: closed %.3f openA %.3f openB %.3f' % (names[k], (OWN[rows] == k).sum(), *T[k]) for k in range(K))), '->', ' '.join('%s=%s' % (r, role[r]) for r in ('c', 'oA', 'oB', 'rest') if r in role))
SA = np.load(a.scores + '/scores_pA.npy', mmap_mode='r')[rows]; SB = np.load(a.scores + '/scores_pB.npy', mmap_mode='r')[rows]
ok = np.isfinite(SA).all(1) & np.isfinite(SB).all(1); rows, SA, SB = rows[ok], np.asarray(SA[ok], np.float64), np.asarray(SB[ok], np.float64)
lam = np.asarray(M['ppnoise'])[rows].astype(np.float64); lam = np.where(np.isfinite(lam) & (lam > 0), lam, np.nanmedian(lam)); SA /= lam[:, None]; SB /= lam[:, None]; own = OWN[rows]
def z(x): return ((x - x.mean()) / x.std())[:, None]
for spec in a.pairs.split(','):
    w = spec.split(':'); tag, ia, ib = w[0], names.index(role[w[1]]), names.index(role[w[2]])
    m = np.isin(own, [names.index(role[s]) for s in w[3].split('+')]) if len(w) > 3 else np.ones(len(rows), bool)
    dA, dB = SA[m, ia] - SA[m, ib], SB[m, ia] - SB[m, ib]; D = os.path.join(a.out, tag, a.name); os.makedirs(D, exist_ok=True)
    np.savez(os.path.join(D, 'coords.npz'), rows=rows[m], tfull=-z(dA + dB), tA=-z(dA), tB=-z(dB), ax=np.array([0]), lam=np.array([np.nan]), nper=np.array([0]))
    open(os.path.join(D, 'arm_final_all.txt'), 'w').write('# orientation stub for lsamp_repro_bins3.py (no maps): pole a is low\nstate_01.mrc 1 0 stub\nstate_02.mrc 0 1 stub\n')
    print('%s %s: %s - %s over %d of %d particles; fold A vs B corr %.3f' % (a.name, tag, role[w[1]], role[w[2]], m.sum(), len(rows), np.corrcoef(dA, dB)[0, 1]), flush=True)
