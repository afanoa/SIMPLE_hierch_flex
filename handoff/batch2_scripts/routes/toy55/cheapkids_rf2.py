#!/usr/bin/env python3
"""cheapkids_rf2.py = cheapkids_rf.py + --axis w1|w2: EQUAL-WIDTH bins of the view-centred PC1 / PC2 over its P..100-P percentile range
(--p 1; the outer bins take the rest): small end bins that sample the extremes. cheapkids_rf.py = cheapkids.py (copied 2026-10-21 from the sub-class session) + --axis 1|2|grid: view-balanced quantile bins along PC1,
PC2, or a GxG grid of the two (G = --k), all inside every direction bin; output file children_<axis>.txt.
Cheap view-balanced per-particle children of one tree node: the node's flex latent of its members (at the parent's poses), centred
per direction bin (spherical k-means of the projection directions, as dir_bins.py), PC1 of the centred latent, K equal quantile bins
of PC1 INSIDE every direction bin (every child takes 1/K of every bin = view-balanced by construction). Writes children.txt (per row,
0 = not a member) and coords.txt (per row PC1 PC2, nan outside). Optional truth: |corr| of PC1/PC2 with a coded motion.
usage: cheapkids.py --oris ptcl3D_dump.txt --rows rows.txt --flex flex_pca_coordinates.txt --out DIR [--k 10] [--perbin 1300]
       [--gt gt --codes swing:200:230]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--oris', required=True); ap.add_argument('--rows', required=True); ap.add_argument('--flex', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--k', type=int, default=10); ap.add_argument('--perbin', type=int, default=1300)
ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='swing:200:230'); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--axis', default='1'); ap.add_argument('--p', type=float, default=1.)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
rows = np.loadtxt(a.rows).astype(int); NR = len(rows); e1 = np.zeros(NR); e2 = np.zeros(NR); i = 0
for ln in open(a.oris):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); e1[i] = float(kv['e1']); e2[i] = float(kv['e2']); i += 1
assert i == NR, 'oris rows %d vs node rows %d' % (i, NR)
F = np.loadtxt(a.flex, comments='#', ndmin=2); fi = F[:, 0].astype(int); keep = rows[fi - 1] > 0; p = fi[keep]; Z = F[keep, 5:]; N = len(p)
E1, E2 = np.radians(e1[p - 1]), np.radians(e2[p - 1]); v = np.c_[np.sin(E2) * np.cos(E1), np.sin(E2) * np.sin(E1), np.cos(E2)]
nb = max(20, int(round(N / a.perbin))); C = v[rng.choice(N, nb, replace=False)].copy()
for it in range(40):   # spherical k-means of the directions (dir_bins.py)
    lab = np.argmax(v @ C.T, 1)
    for k in range(nb):
        m = lab == k
        C[k] = v[m].mean(0) / np.linalg.norm(v[m].mean(0)) if m.any() else v[rng.integers(N)]
lab = np.argmax(v @ C.T, 1); _, bi = np.unique(lab, return_inverse=True); nb = bi.max() + 1; cnt = np.bincount(bi)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=nb) for j in range(Z.shape[1])], 1) / cnt[:, None]; Zc = Z - bm[bi]; Zc -= Zc.mean(0)
w, U = np.linalg.eigh(np.cov(Zc.T)); U = U[:, ::-1]; w = w[::-1]; y = Zc @ U[:, 0]; y2 = Zc @ U[:, 1]
ch = np.zeros(N, int)
def qbins(x, m, k): r = np.argsort(np.argsort(x[m])); return np.minimum(r * k // len(m), k - 1)
if a.axis in ('w1', 'w2'):
    x = y if a.axis == 'w1' else y2; lo, hi = np.percentile(x, [a.p, 100 - a.p]); e = np.linspace(lo, hi, a.k + 1)
    ch = np.clip(np.searchsorted(e, x, side='right') - 1, 0, a.k - 1) + 1
else:
  for b in range(nb):
    m = np.where(bi == b)[0]
    if a.axis == 'grid': ch[m] = qbins(y, m, a.k) * a.k + qbins(y2, m, a.k) + 1
    else: ch[m] = qbins(y if a.axis == '1' else y2, m, a.k) + 1
out = np.zeros(NR, int); out[p - 1] = ch; np.savetxt(os.path.join(a.out, 'children_%s.txt' % a.axis), out, fmt='%d')
co = np.full((NR, 2), np.nan); co[p - 1, 0] = y; co[p - 1, 1] = y2; np.savetxt(os.path.join(a.out, 'coords.txt'), co, fmt='%.5f')
print('node: %d members, %d direction bins (members min %d median %d), %d-dim latent; PC variance shares of the bin-centred latent %s'
      % (N, nb, cnt.min(), np.median(cnt), Z.shape[1], ' '.join('%.2f' % s for s in (w / w.sum())[:5])))
print('children (axis %s): populations %s' % (a.axis, ' '.join(str((ch == c).sum()) for c in range(1, ch.max() + 1))))
if a.gt:
    G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(NR + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
    nm, lo, hi = a.codes.split(':'); s = np.where((conf >= int(lo)) & (conf <= int(hi)), conf - int(lo), np.nan); d = np.isfinite(s)
    if d.sum() > 100: print('TRUTH %s particles %d (%.0f %%); |corr| per particle with PC1 %.3f, PC2 %.3f' % (nm, d.sum(), 100 * d.mean(), abs(np.corrcoef(y[d], s[d])[0, 1]), abs(np.corrcoef(y2[d], s[d])[0, 1])))
