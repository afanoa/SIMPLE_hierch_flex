#!/usr/bin/env python3
"""State graph of a final labelling (both halves, one label per project row, 0 = none): nodes = states (size = particles,
colour = truth majority when GT is given), edge weight = cross-half similarity of the state deviations from the
consensus (39-10 A band: state i in one half against state j in the other, averaged both ways), dashed edges = particle
sharing (cosine overlap of the soft memberships from one cross-half E-step, consensus + FSC-weighted deviations, <= 20 A).
Layout = classical MDS of 1 - similarity; right panel = similarity matrix ordered by average linkage with its dendrogram.
usage: state_graph.py --sim A,B --poses A,B --labels A,B --cache PREFIX --png OUT [--gt A,B --truth NPZ --ctfabs 1]"""
import sys, os, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
ap = argparse.ArgumentParser()
for k in ('sim', 'poses', 'labels', 'cache', 'png'): ap.add_argument('--' + k, required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--cut', type=float, default=0.45); ap.add_argument('--title', default='')
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7)
ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--lim', type=float, default=20.); ap.add_argument('--bandlo', type=float, default=40.)
ap.add_argument('--ctfabs', type=int, default=0, help='1: |CTF| (phase-flipped stacks)'); ap.add_argument('--minsize', type=int, default=200); a = ap.parse_args()
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def load_half(simdir, posef, labf, tag):   # splitbench3.py geometry with its calibrated conv 1,1,1,1,1,1 (transpose: R k)
    lab = np.loadtxt(labf).astype(int); rows = np.where(lab > 0)[0]
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    n = int(np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    F = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r')[rows]; pp = P[rows]; kf = np.hypot(fx2, fy2) / a.smpd; ang = np.arctan2(fy2, fx2)
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
    if a.ctfabs: C = np.abs(C)
    ramp = np.exp(2j * np.pi * (fx2[None] * pp[:, 3, None, None] + fy2[None] * pp[:, 4, None, None])).astype(np.complex64)
    Wv = (C * np.conj(ramp) * F).reshape(len(rows), -1); Wd = (C ** 2).reshape(len(rows), -1); F2 = (np.abs(F) ** 2).reshape(len(rows), -1).astype(np.float32)
    q = np.round(np.einsum('qj,nij->nqi', kk0, rotm(pp[:, 0], pp[:, 1], pp[:, 2]))).astype(np.int32); ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
    idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3
    log('%s: %d particles in states %s' % (tag, len(rows), np.bincount(lab[rows]).tolist()[1:])); return dict(lab=lab[rows], rows=rows, Wv=Wv, Wd=Wd, F2=F2, idx=idx, n=n, kk0=kk0)
HA = load_half(*[s.split(',')[0] for s in (a.sim, a.poses, a.labels)], 'A'); HB = load_half(*[s.split(',')[1] for s in (a.sim, a.poses, a.labels)], 'B')
c = a.crop; n = HA['n']; m3 = c ** 3 + 1; ns = c // 2 + 1; khi = c // 2 - 1; klo = max(int(round(n * a.smpd / a.bandlo)), 3)
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel()
g_ = np.arange(c) - c // 2; r3 = np.sqrt(g_[:, None, None] ** 2 + g_[None, :, None] ** 2 + g_[None, None, :] ** 2) * n / c; msk3 = np.clip((230. / 2 / a.smpd + 4 - r3) / 8., 0, 1)
ids = sorted(set(HA['lab']) & set(HB['lab'])); ids = [k for k in ids if (HA['lab'] == k).sum() >= a.minsize and (HB['lab'] == k).sum() >= a.minsize]; K = len(ids)
for H in (HA, HB): H['own'] = np.array([ids.index(v) if v in ids else -1 for v in H['lab']])
def acc(H, sel):
    ii = H['idx'][sel].ravel(); return np.bincount(ii, H['Wv'][sel].real.ravel(), m3) + 1j * np.bincount(ii, H['Wv'][sel].imag.ravel(), m3), np.bincount(ii, H['Wd'][sel].ravel(), m3)
ACC = {h: [acc(H, np.where(H['own'] == j)[0]) + (int((H['own'] == j).sum()),) for j in range(K)] for h, H in ((0, HA), (1, HB))}
nref = int(np.median([v[2] for v in ACC[0]])); d0 = ACC[0][0][1][:-1] * nref / ACC[0][0][2]
epsr = np.array([np.median(d0[(sh == k) & (d0 > 0)]) if ((sh == k) & (d0 > 0)).any() else 1. for k in range(ns)]) * a.eps; epsv = np.r_[epsr[sh], epsr[-1]]
def vmap(h, js, mask=True):
    num = sum(ACC[h][j][0] for j in js); den = sum(ACC[h][j][1] for j in js); N = sum(ACC[h][j][2] for j in js); V = num / (den + epsv * N / nref)
    if not mask: return V
    v = np.fft.fftshift(np.real(np.fft.ifftn(V[:-1].reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
# ---- (1) similarity of the state deviations across halves
cons = [vmap(h, range(K)) for h in (0, 1)]; D = [[vmap(h, [j]) - cons[h] for j in range(K)] for h in (0, 1)]
band = (sh >= klo) & (sh <= khi); coarse = (sh >= 2) & (sh < klo)
def corr(X, Y, m): x, y = X[m], Y[m]; return float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30))
S = np.array([[0.5 * (corr(D[0][i], D[1][j], band) + corr(D[1][i], D[0][j], band)) for j in range(K)] for i in range(K)])
Sc = np.array([[0.5 * (corr(D[0][i], D[1][j], coarse) + corr(D[1][i], D[0][j], coarse)) for j in range(K)] for i in range(K)])
# ---- (2) particle sharing: one cross-half E-step (consensus + deviation x cross-half FSC per shell, <= lim A)
kb = min(c // 2 - 1, int(round(n * a.smpd / a.lim))); r2 = np.round(np.hypot(HA['kk0'][:, 0], HA['kk0'][:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]
U = [[vmap(h, [j], mask=False) for j in range(K)] for h in (0, 1)]; UC = [vmap(h, range(K), mask=False) for h in (0, 1)]; REF = {0: [], 1: []}
for j in range(K):
    DA, DB = U[0][j][:-1] - UC[0][:-1], U[1][j][:-1] - UC[1][:-1]
    f = np.bincount(sh, np.real(DA * np.conj(DB)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(DA) ** 2, ns) * np.bincount(sh, np.abs(DB) ** 2, ns)), 1e-30)
    f = np.clip(np.convolve(np.pad(f, 1, mode='edge'), np.ones(3) / 3., mode='valid'), 0., 1.)
    REF[0].append(np.r_[UC[0][:-1] + f[sh] * DA, 0.].astype(np.complex64)); REF[1].append(np.r_[UC[1][:-1] + f[sh] * DB, 0.].astype(np.complex64))
def estep(H, R):
    ix = H['idx'][:, cols]; Wv = H['Wv'][:, cols]; Wd = H['Wd'][:, cols]; F2 = H['F2'][:, cols]; shq = np.minimum(sh[ix % (c ** 3)], ns - 1); own = H['own']; R2 = np.zeros(ns); C2 = np.zeros(ns)
    for j in range(K):
        r = np.where(own == j)[0]; Pg = R[j][ix[r]]; res = F2[r] - 2. * np.real(Pg * np.conj(Wv[r])) + Wd[r] * np.abs(Pg) ** 2; R2 += np.bincount(shq[r].ravel(), res.ravel(), ns); C2 += np.bincount(shq[r].ravel(), minlength=ns)
    winv = (1. / np.maximum(R2 / np.maximum(C2, 1), 1e-12)).astype(np.float32)[shq]; Sc_ = np.zeros((len(own), K), np.float32)
    for j in range(K):
        for b in range(0, len(own), 8000):
            s_ = slice(b, b + 8000); Pg = R[j][ix[s_]]; Sc_[s_, j] = 0.5 * ((2. * np.real(Pg * np.conj(Wv[s_])) - Wd[s_] * np.abs(Pg) ** 2) * winv[s_]).sum(1)
    pri = np.array([(own == j).mean() for j in range(K)]); Sc_ += np.log(np.maximum(pri, 1e-6))[None]; Sc_ -= Sc_.max(1, keepdims=True); W = np.exp(Sc_); return W / W.sum(1, keepdims=True)
O = np.zeros((K, K)); WS = []
for H, R in ((HA, REF[1]), (HB, REF[0])):
    W = estep(H, R); G = W.T @ W; O += 0.5 * G / np.sqrt(np.outer(np.diag(G), np.diag(G))); WS.append(W.astype(np.float16))
# ---- truth annotation
ann = ['' for _ in range(K)]; col = ['#999999'] * K
if a.gt:
    cf = [np.loadtxt(f, comments='#')[:, 1].astype(int) for f in a.gt.split(',')]; tops = []
    clade = None
    if a.truth: z = np.load(a.truth); clade = dict(zip(z['ids'], fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=a.cut, criterion='distance')))
    for j in range(K):
        t = np.r_[cf[0][HA['rows'][HA['own'] == j]], cf[1][HB['rows'][HB['own'] == j]]]; real = t[t > 0]
        if (t == 0).mean() > 0.5: ann[j] = 'junk %.0f%%' % (100 * (t == 0).mean()); tops.append(0); continue
        bc = np.bincount(real); tp = int(bc.argmax()); tops.append(tp); ann[j] = 'c%d %.0f%%' % (tp, 100. * bc[tp] / len(real))
    pal = plt.cm.tab10(np.arange(10)); u = sorted(set(tops)); cm = {v: ('#bbbbbb' if v == 0 else pal[i % 10]) for i, v in enumerate(u)}; col = [cm[v] for v in tops]
size = np.array([ACC[0][j][2] + ACC[1][j][2] for j in range(K)], float)
log('states %s; sizes %s' % (ids, size.astype(int).tolist()))
log('cross-half self-reproducibility (diagonal, 39-10 A): %s' % ' '.join('%.2f' % v for v in np.diag(S)))
for i in range(K):
    nb = [j for j in np.argsort(-S[i]) if j != i][:2]; ob = [j for j in np.argsort(-O[i]) if j != i][:2]
    log('  state %d (%s, n %d): most similar %s; most shared particles %s' % (ids[i], ann[i], size[i], ', '.join('%d (%.2f)' % (ids[j], S[i, j]) for j in nb), ', '.join('%d (%.2f)' % (ids[j], O[i, j]) for j in ob)))
# ---- figure
Dm = 1. - 0.5 * (S + S.T); np.fill_diagonal(Dm, 0.); Dm = np.maximum(Dm, 0.)
J = np.eye(K) - 1. / K; B = -0.5 * J @ (Dm ** 2) @ J; ev, evec = np.linalg.eigh(B); o = np.argsort(ev)[::-1][:2]; XY = evec[:, o] * np.sqrt(np.maximum(ev[o], 1e-9))
fig = plt.figure(figsize=(15, 6.8)); ax = fig.add_axes([0.02, 0.06, 0.5, 0.86]); ax2 = fig.add_axes([0.6, 0.06, 0.36, 0.62]); axd = fig.add_axes([0.6, 0.7, 0.36, 0.22])
offd = S[~np.eye(K, dtype=bool)]; thr = max(0.1, np.quantile(offd, 0.6)) if len(offd) else 0.1
for i in range(K):
    for j in range(i + 1, K):
        s = 0.5 * (S[i, j] + S[j, i])
        if s > thr: ax.plot(XY[[i, j], 0], XY[[i, j], 1], '-', color='#3b6ea5', lw=1 + 6 * (s - thr) / max(1 - thr, 1e-6), alpha=0.7, zorder=1)
        if O[i, j] > 0.15: ax.plot(XY[[i, j], 0], XY[[i, j], 1], '--', color='#c0504d', lw=1 + 5 * O[i, j], alpha=0.8, zorder=1)
ax.scatter(XY[:, 0], XY[:, 1], s=300 + 2500 * size / size.max(), c=col, edgecolors='k', zorder=2)
for i in range(K): ax.text(XY[i, 0], XY[i, 1], '%d\n%s' % (ids[i], ann[i]), ha='center', va='center', fontsize=8, zorder=3)
ax.set_title('state graph: blue = cross-half map similarity (> %.2f), red dashed = shared particles (> 0.15); area = particles' % thr, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
pad = 0.25 * (np.ptp(XY, 0).max() + 1e-3); ax.set_xlim(XY[:, 0].min() - pad, XY[:, 0].max() + pad); ax.set_ylim(XY[:, 1].min() - pad, XY[:, 1].max() + pad)
Z = linkage(squareform(Dm, checks=False), 'average'); dn = dendrogram(Z, ax=axd, no_labels=True, color_threshold=0, above_threshold_color='k'); order = dn['leaves']; axd.set_yticks([]); axd.set_title('average linkage of 1 - similarity (drawn after the fact, not a sorting order)', fontsize=9)
for sp in axd.spines.values(): sp.set_visible(False)
v = np.abs(S).max(); im = ax2.imshow(S[np.ix_(order, order)], cmap='RdBu_r', vmin=-v, vmax=v)
ax2.set_xticks(range(K)); ax2.set_xticklabels([str(ids[i]) for i in order], fontsize=8); ax2.set_yticks(range(K)); ax2.set_yticklabels(['%d %s' % (ids[i], ann[i]) for i in order], fontsize=8)
plt.colorbar(im, ax=ax2, fraction=0.04); fig.suptitle(a.title or 'state graph of the final states', fontsize=11)
plt.savefig(a.png, dpi=100, bbox_inches='tight')
cfs = [np.loadtxt(f, comments='#')[:, 1].astype(int)[H['rows']] for f, H in zip(a.gt.split(','), (HA, HB))] if a.gt else [np.zeros(len(HA['rows']), int), np.zeros(len(HB['rows']), int)]
np.savez(a.png.replace('.png', '.npz'), ids=ids, S=S, Sc=Sc, O=O, size=size, xy=XY, WA=WS[0], WB=WS[1], ownA=HA['own'], ownB=HB['own'], confA=cfs[0], confB=cfs[1]); log('png', a.png)
