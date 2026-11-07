#!/usr/bin/env python3
"""lsamp_pq_coords.py: pole-pair coordinates of one run for a node with three poles (closed c, open A, open B), from a score space that
holds the three pole maps (per pixel set, noise-normalised as lsamp2h.py --ppnorm 1): c_oA = s(c) - s(oA), c_oB = s(c) - s(oB),
oA_oB = s(oA) - s(oB), c_rest = s(c) - log(exp s(oA) + exp s(oB)) (closed against the other two). Each is standardised and written as
OUT/<type>/<name>/coords.npz (tA, tB, tfull; first-named pole low) with an orientation stub for lsamp_repro_bins2.py. --roles auto names the
poles of a score space by its own leaves' overlap with reference posteriors (--ref_post closed,openA,openB: master-row posteriors).
usage: lsamp_pq_coords.py --scores SC --name RUN --out OUT --roles C,OA,OB | auto --ref_post P1,P2,P3 [--types c_oA,c_oB,oA_oB,c_rest]"""
import os, argparse, itertools, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--name', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--roles', required=True); ap.add_argument('--ref_post', default=''); ap.add_argument('--types', default='c_oA,c_oB,oA_oB,c_rest')
a = ap.parse_args()
M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; rows = np.where((OWN >= 0) & (EO >= 0))[0]
if a.roles == 'auto':   # name the poles: the assignment of leaves to (closed, open A, open B) with the largest summed mean reference posterior
    R = [np.loadtxt(f) for f in a.ref_post.split(',')]; K = len(names); T = np.array([[R[j][rows[OWN[rows] == k]].mean() for j in range(3)] for k in range(K)])
    best = max(itertools.permutations(range(K), 3), key=lambda p: sum(T[p[j], j] for j in range(3))); roles = [names[k] for k in best]
    print('%s roles: mean reference posterior of each leaf\'s particles (closed / open A / open B): %s -> c %s, oA %s, oB %s' % (a.name, '; '.join('%s %d: %s' % (names[k], (OWN[rows] == k).sum(), ' '.join('%.3f' % v for v in T[k])) for k in range(K)), *roles))
else: roles = a.roles.split(',')
ic, ia, ib = [names.index(r) for r in roles]
SA = np.load(a.scores + '/scores_pA.npy', mmap_mode='r')[rows]; SB = np.load(a.scores + '/scores_pB.npy', mmap_mode='r')[rows]
ok = np.isfinite(SA).all(1) & np.isfinite(SB).all(1); rows, SA, SB = rows[ok], np.asarray(SA[ok], np.float64), np.asarray(SB[ok], np.float64)
lam = np.asarray(M['ppnoise'])[rows].astype(np.float64); lam = np.where(np.isfinite(lam) & (lam > 0), lam, np.nanmedian(lam)); SA /= lam[:, None]; SB /= lam[:, None]
def coord(S, t):
    if t == 'c_oA': return S[:, ic] - S[:, ia]
    if t == 'c_oB': return S[:, ic] - S[:, ib]
    if t == 'oA_oB': return S[:, ia] - S[:, ib]
    return S[:, ic] - np.logaddexp(S[:, ia], S[:, ib])
def z(x): return ((x - x.mean()) / x.std())[:, None]
for t in a.types.split(','):
    D = os.path.join(a.out, t, a.name); os.makedirs(D, exist_ok=True); dA, dB, dF = coord(SA, t), coord(SB, t), coord(SA + SB, t)
    np.savez(os.path.join(D, 'coords.npz'), rows=rows, tfull=-z(dF), tA=-z(dA), tB=-z(dB), ax=np.array([0]), lam=np.array([np.nan]), nper=np.array([0]))
    open(os.path.join(D, 'arm_final_all.txt'), 'w').write('# orientation stub for lsamp_repro_bins2.py (no maps): the first-named pole is low\nstate_01.mrc 1 0 stub\nstate_02.mrc 0 1 stub\n')
    print('%s %s: %d particles, poles %s; fold A vs B corr %.3f' % (a.name, t, len(rows), '/'.join(roles), np.corrcoef(dA, dB)[0, 1]), flush=True)
