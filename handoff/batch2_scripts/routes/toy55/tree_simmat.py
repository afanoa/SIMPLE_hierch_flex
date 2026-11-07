#!/usr/bin/env python3
"""Leaf x leaf similarity of a data-derived tree, ordered by the tree: entry (i, j) is the cross-half correlation of leaf i's
deviation from the consensus in one half with leaf j's in the other (averaged both ways), on the coarse (>bandlo) and the tested
band. Blocks along the diagonal = branches whose leaves resemble each other; strong entries off the blocks = leaves that sit
closer to another branch. Prints, per leaf, its most similar other leaf and whether it is in the same first-level branch.
usage: tree_simmat.py --sim A,B --poses A,B --leaves A,B --cache PREFIX --png OUT [--gt A,B --truth NPZ]"""
import sys, os, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser()
ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--leaves', required=True); ap.add_argument('--cache', required=True); ap.add_argument('--png', required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--cut', type=float, default=0.45); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--title', default=''); ap.add_argument('--exclude', default='', help='branches (comma list, e.g. r.3) left out of the reference map, e.g. the junk branch')
ap.add_argument('--bandlo', type=float, default=40.); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--eps', type=float, default=0.05)
a = ap.parse_args()
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def half_acc(simdir, posef, leaff, tag):   # per-leaf Fourier accumulators (splitbench3.py geometry, identity conventions)
    lab = np.loadtxt(leaff).astype(int); names = {int(l.split()[0]): l.split()[1] for l in open(leaff + '.names') if l.strip()}
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    nx = int(np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop; n = nx
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    Fc = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r'); kf = np.hypot(fx2, fy2) / a.smpd; ang = np.arctan2(fy2, fx2); out = {}
    for g, nm in names.items():
        rows = np.where(lab == g)[0]
        if not len(rows): continue
        num = np.zeros(c ** 3 + 1, np.complex128); den = np.zeros(c ** 3 + 1)
        for b in range(0, len(rows), 4000):
            r = rows[b:b + 4000]; pp = P[r]; F = np.asarray(Fc[r])
            df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
            chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))
            ramp = np.exp(2j * np.pi * (fx2[None] * pp[:, 3, None, None] + fy2[None] * pp[:, 4, None, None])); Wv = (C * np.conj(ramp) * F).reshape(len(r), -1); Wd = (C ** 2).reshape(len(r), -1)
            q = np.round(np.einsum('qj,nij->nqi', kk0, rotm(pp[:, 0], pp[:, 1], pp[:, 2]))).astype(np.int64); ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
            ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~ok] = c ** 3; ii = ii.ravel()
            num += np.bincount(ii, Wv.real.ravel(), c ** 3 + 1) + 1j * np.bincount(ii, Wv.imag.ravel(), c ** 3 + 1); den += np.bincount(ii, Wd.ravel(), c ** 3 + 1)
        out[nm] = (num[:-1], den[:-1], len(rows))
    log('%s: %d leaves, %d particles' % (tag, len(out), sum(v[2] for v in out.values()))); return out, lab, names
AA, labA, nmA = half_acc(*[s.split(',')[0] for s in (a.sim, a.poses, a.leaves)], 'A'); BB, labB, nmB = half_acc(*[s.split(',')[1] for s in (a.sim, a.poses, a.leaves)], 'B')
nx = int(np.frombuffer(open(a.sim.split(',')[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop; n = nx
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel(); ns = c // 2 + 1
g = np.arange(c) - c // 2; r3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) * n / c; msk3 = np.clip((a.mskdiam / 2 / a.smpd + 4 - r3) / 8., 0, 1)
klo = max(int(round(n * a.smpd / a.bandlo)), 3); khi = c // 2 - 1
key = lambda s: [int(x) for x in s.split('.')[1:]]; leaves = sorted([L for L in AA if L in BB], key=key)
nref = int(np.median([AA[L][2] for L in leaves])); den0 = 0.5 * (AA[leaves[0]][1] + BB[leaves[0]][1]) * nref / AA[leaves[0]][2]
epsr = np.array([np.median(den0[(sh == k) & (den0 > 0)]) if ((sh == k) & (den0 > 0)).any() else 1. for k in range(ns)]) * a.eps
def vmap(D, Ls):
    num = sum(D[L][0] for L in Ls); den = sum(D[L][1] for L in Ls); N = sum(D[L][2] for L in Ls); V = num / (den + epsr[sh] * N / nref)
    v = np.fft.fftshift(np.real(np.fft.ifftn(V.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
xl = [x for x in a.exclude.split(',') if x]; refl = [L for L in leaves if not any(L == x or L.startswith(x + '.') for x in xl)]; cA, cB = vmap(AA, refl), vmap(BB, refl); log('reference map from %d of %d leaves%s' % (len(refl), len(leaves), ' (excluded %s)' % a.exclude if xl else '')); MA = {L: vmap(AA, [L]) - cA for L in leaves}; MB = {L: vmap(BB, [L]) - cB for L in leaves}
bands = {'coarse >%.0f A' % a.bandlo: (sh >= 2) & (sh < klo), 'tested %.0f-%.1f A' % (n * a.smpd / klo, n * a.smpd / khi): (sh >= klo) & (sh <= khi)}
def corr(X, Y, m): x, y = X[m], Y[m]; return float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30))
S = {b: np.array([[0.5 * (corr(MA[i], MB[j], m) + corr(MB[i], MA[j], m)) for j in leaves] for i in leaves]) for b, m in bands.items()}
first = np.array(['.'.join(L.split('.')[:2]) for L in leaves]); ann = {L: '' for L in leaves}
if a.gt:   # truth annotation: majority conformer (and reference clade) and junk share, both halves pooled
    clade = None
    if a.truth:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        z = np.load(a.truth); clade = dict(zip(z['ids'], fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=a.cut, criterion='distance')))
    cf = [np.loadtxt(f, comments='#')[:, 1].astype(int) for f in a.gt.split(',')]
    for L in leaves:
        t = np.r_[cf[0][labA == {v: k for k, v in nmA.items()}[L]], cf[1][labB == {v: k for k, v in nmB.items()}[L]]]; real = t[t > 0]
        if not len(real): ann[L] = ' junk'; continue
        bc = np.bincount(real); top = int(bc.argmax()); s = ' c%d %.0f%%' % (top, 100. * bc[top] / len(real))
        if clade is not None: cc = np.bincount([clade.get(k, 0) for k in real]); s += ' clade%d %.0f%%' % (cc.argmax(), 100. * cc.max() / len(real))
        ann[L] = s + (' junk %.0f%%' % (100. * (t == 0).mean()) if (t == 0).mean() > 0.1 else '')
for b, M in S.items():
    off = M - np.diag(np.diag(M)) - 9 * np.eye(len(leaves)); nn = off.argmax(1); same = first[nn] == first
    if xl:   # nearest-leaf readout among the leaves of the reference only
        ex = np.array([L not in refl for L in leaves]); off[:, ex] = -9
    wb = M[(first[:, None] == first[None]) & ~np.eye(len(leaves), dtype=bool)]; bb = M[first[:, None] != first[None]]
    log('\n[%s] diagonal (own cross-half reproducibility) median %.2f; within-branch mean %.2f, between-branch mean %.2f; nearest other leaf in the same first-level branch for %d/%d leaves' % (b, np.median(np.diag(M)), wb.mean() if len(wb) else np.nan, bb.mean(), same.sum(), len(leaves)))
    for i, L in enumerate(leaves): log('  %-10s self %.2f  nearest %-10s %.2f %s%s' % (L, M[i, i], leaves[nn[i]], M[i, nn[i]], '' if same[i] else 'OTHER BRANCH', ann[L]))
fig, axs = plt.subplots(1, 2, figsize=(2 * (0.36 * len(leaves) + 3.2), 0.36 * len(leaves) + 2.2))
for ax, (b, M) in zip(axs, S.items()):
    v = np.abs(M).max(); im = ax.imshow(M, cmap='RdBu_r', vmin=-v, vmax=v); ax.set_title(b, fontsize=10)
    ax.set_xticks(range(len(leaves))); ax.set_xticklabels(leaves, rotation=90, fontsize=7); ax.set_yticks(range(len(leaves))); ax.set_yticklabels(['%s%s' % (L, ann[L]) for L in leaves], fontsize=7)
    for k in range(1, len(leaves)):
        if first[k] != first[k - 1]: ax.axhline(k - 0.5, color='k', lw=1); ax.axvline(k - 0.5, color='k', lw=1)
    plt.colorbar(im, ax=ax, fraction=0.04)
fig.suptitle(a.title or 'leaf x leaf cross-half similarity (deviation from the %s), tree order; lines = first-level branches' % ('mean of the leaves outside ' + a.exclude if xl else 'consensus'), fontsize=10)
plt.tight_layout(); plt.savefig(a.png, dpi=100, bbox_inches='tight'); np.savez(a.png.replace('.png', '.npz'), leaves=leaves, **{('S%d' % i): M for i, M in enumerate(S.values())}); log('png', a.png)
