#!/usr/bin/env python3
"""Sampling the final states (after the merge-up): every final state is a group of leaves; its half-maps are the particle-weighted
averages of its leaves' even / odd maps. The states are organised by a tree (complete linkage on the cross-half relative difference,
40-10 A) and, at every tree node holding >= MINST states (the root first), the principal axes of the node's state maps come from the
particle-weighted, centred cross-half Gram matrix (unbiased); the leading axes up to SHARE of the node's variance (at most MAXAX) are
reported with the node's states in order along each axis and its two ends. With GT + CODES: per axis the Spearman correlation of the
coordinate with the per-state mean of every truth coordinate, and the truth at the two ends.
usage: state_axes_tree.py OUT --rdir LEAF_HALVES_DIR --leaves LEAF_LABELS --units FINAL_UNITS --unit_leaves FINAL_LEAVES
       [--unfil --hp 40 --lp 10 --minst 3 --maxax 3 --share 0.9 --gt GT --codes swing:200:1]"""
import os, sys, argparse, numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', required=True); ap.add_argument('--leaves', required=True); ap.add_argument('--units', required=True)
ap.add_argument('--unit_leaves', required=True); ap.add_argument('--unfil', action='store_true'); ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lp', type=float, default=10.)
ap.add_argument('--minst', type=int, default=3); ap.add_argument('--maxax', type=int, default=3); ap.add_argument('--share', type=float, default=0.9); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
lab = np.loadtxt(a.leaves, dtype=int); units = np.loadtxt(a.units, dtype=int); sfx = '_unfil' if a.unfil else ''
UL = {int(l.split()[0]): [int(x) for x in l.split()[1].split(',')] for l in open(a.unit_leaves) if l.strip() and not l.startswith('#')}
U = sorted(u for u in UL if (units == u).sum() > 0); nl = {lf: float((lab == lf).sum()) for u in U for lf in UL[u]}
p0 = '%s/recvol_state%02d_even%s.mrc' % (a.rdir, UL[U[0]][0], sfx); n = rd(p0).shape[0]; smpd = float(np.fromfile(p0, np.float32, 13)[10]) / n
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9); band = (res <= a.hp) & (res >= a.lp); wgt = np.where((kz[band] == 0) | (np.abs(kz[band]) == 0.5), 1., 2.)
x = np.arange(n) - n / 2.; M = (np.sqrt(sum(g ** 2 for g in np.meshgrid(x, x, x, indexing='ij'))) < 0.45 * n).astype(float)   # spherical support (no data-specific mask)
E, O, w = [], [], []
for u in U:
    ws = np.array([nl[lf] for lf in UL[u]]); ws /= ws.sum()
    E.append(sum(wi * np.fft.rfftn(M * rd('%s/recvol_state%02d_even%s.mrc' % (a.rdir, lf, sfx)))[band] for wi, lf in zip(ws, UL[u])))
    O.append(sum(wi * np.fft.rfftn(M * rd('%s/recvol_state%02d_odd%s.mrc' % (a.rdir, lf, sfx)))[band] for wi, lf in zip(ws, UL[u]))); w.append((units == u).sum())
E, O, w = np.array(E), np.array(O), np.array(w, float); C = np.real((E * wgt) @ O.conj().T); C = (C + C.T) / 2; N = len(U)
D = np.zeros((N, N))
for i in range(N):
    for j in range(i + 1, N): D[i, j] = D[j, i] = (C[i, i] + C[j, j] - C[i, j] - C[j, i]) / np.sqrt(max(C[i, i], 1e-30) * max(C[j, j], 1e-30))
iu = np.triu_indices(N, 1); Z = linkage(np.maximum(D[iu], 0), method='complete'); root = to_tree(Z)
truth = {}
if a.gt:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(units), -1); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
    for c in [w_ for w_ in a.codes.split(',') if w_]:
        p = c.split(':')
        if p[0] == 'grid':
            i0, ks, js = int(p[1]), float(p[2]), float(p[3]); m = conf >= i0; truth['grid_k'] = np.where(m, ((conf - i0) // 10) * ks, np.nan); truth['grid_j'] = np.where(m, ((conf - i0) % 10) * js, np.nan)
        else:
            nm, i0, st = p[0], int(p[1]), float(p[2]); hi = min([int(q.split(':')[1]) for q in a.codes.split(',') if q and not q.startswith('grid:') and int(q.split(':')[1]) > i0] + [10 ** 9])
            truth[nm] = np.where((conf >= i0) & (conf < hi), (conf - i0) * st, np.nan)
def tmean(u, nm):
    v = truth[nm][units == u]; v = v[np.isfinite(v)]; return (v.mean(), len(v) / max((units == u).sum(), 1)) if len(v) else (np.nan, 0.)
rk = lambda v: np.argsort(np.argsort(v)).astype(float)
def axes(node, path):
    ids = node.pre_order(); sts = [U[i] for i in ids]
    if len(ids) < a.minst: return
    wn = w[ids] / w[ids].sum(); Cn = C[np.ix_(ids, ids)]; m_ = Cn @ wn; Gc = Cn - m_[:, None] - m_[None, :] + wn @ Cn @ wn   # centred Gram
    sw = np.sqrt(wn); lam, V = np.linalg.eigh(sw[:, None] * Gc * sw[None, :]); lam, V = lam[::-1], V[:, ::-1]; pos = lam[lam > 0].sum()
    if pos <= 0: return
    X = V * np.sqrt(np.maximum(lam, 0))[None, :] / sw[:, None]; cum = np.cumsum(lam) / pos; nax = min(a.maxax, int(np.searchsorted(cum, a.share) + 1), len(ids) - 1)
    print('== node %s: %d states (%d particles); axis variance shares %s -> %d axes' % (path or 'root', len(ids), w[ids].sum(), ' '.join('%.2f' % (l / pos) for l in lam[:5]), nax))
    for k in range(nax):
        o = np.argsort(X[:, k]); line = '   axis %d (%.0f %%): %s' % (k + 1, 100 * lam[k] / pos, ' '.join('%d' % sts[i] for i in o))
        for nm in truth:
            tv = np.array([tmean(s, nm) for s in sts]); ok = (tv[:, 1] >= 0.5)
            if ok.sum() >= 3: line += ' | %s rho %+.2f, ends %.1f / %.1f' % (nm, np.corrcoef(rk(X[ok, k]), rk(tv[ok, 0]))[0, 1], tv[o[0], 0] if tv[o[0], 1] >= 0.5 else np.nan, tv[o[-1], 0] if tv[o[-1], 1] >= 0.5 else np.nan)
        print(line)
    if node.left: axes(node.left, path + 'L')
    if node.right: axes(node.right, path + 'R')
print('%d final states, %d leaves; tree over the states (complete linkage, relative cross-half difference)' % (N, len(nl)))
axes(root, '')
