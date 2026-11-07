#!/usr/bin/env python3
"""Focus masks for every node of a state tree, from the data. For node u with member leaves j (weights = particles):
    v_u(x) = sum_j w_j (E_j - E_u)(O_j - O_u) / sum_j w_j       (low-passed half-maps, E_u / O_u the members' weighted means)
is the cross-half variance map of the node: unbiased (the halves' noise is independent), zero in expectation wherever the members
agree, and equal to sum_k lambda_k u_k^E(x) u_k^O(x) over the node's principal components. Where nothing varies it is symmetric
about zero, so its negative side measures the noise: the mask keeps voxels above the (1 - ALPHA) quantile of -v_u inside the
sphere, dilated 3 px, soft 2 px. Per node also the cross-half PCA of the members (reproducible components: eigenvalues above the
magnitude of the most negative one) and the overlap of the two leading components' supports (separate supports = two independent
changes in one node). With --gt/--truth: the same variance map of the members' truth maps; the mask's recall (share of the true
variance inside it) and its size against the true support (voxels holding 90 % of the true variance).
node_masks2.py: before the variance, every member's half-maps are regressed per Fourier shell on the node's mean (cross-half
coefficient b_j(s) = Re<E_j, O_u>_s / Re<E_u, O_u>_s; scale, resolution and junk-blur differences removed); masks by three rules:
negative tail (ALPHA), 80 % of the positive variance mass, 10 % of the maximum; each with its recall of the true variance.
usage: node_masks2.py OUT --rdir RD --labels refined_states.txt --tree merge_tree.txt [--gt GT --truth dir1,dir2] [--mskdiam 220
       --lpmask 20 --alpha 0.001 --minleaves 2]"""
import os, glob, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', required=True); ap.add_argument('--labels', required=True)
ap.add_argument('--tree', required=True); ap.add_argument('--gt', default=''); ap.add_argument('--truth', default='')
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--lpmask', type=float, default=20.); ap.add_argument('--alpha', type=float, default=0.001)
ap.add_argument('--minleaves', type=int, default=2)
a = ap.parse_args(); os.makedirs(a.out + '/masks', exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
lab = np.loadtxt(a.labels, dtype=int)
hdr = open(a.tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; NL = len(names)
S = np.atleast_2d(np.loadtxt(a.tree, comments='#')); mem = {i: [i] for i in range(NL)}; par = {}; kids = {}
for st in S: p, q, r = int(st[1]), int(st[2]), int(st[3]); mem[r] = mem[p] + mem[q]; par[p] = par[q] = r; kids[r] = (p, q)
root = int(S[-1, 3])
f0 = '%s/recvol_state%02d_even_unfil.mrc' % (a.rdir, names[0]); n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
lpf = np.exp(-(kr / smpd * a.lpmask) ** 2 * np.log(2)); lp = lambda x: np.fft.irfftn(np.fft.rfftn(x) * lpf, s=(n, n, n))
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
gk, ball = np.exp(-Rw ** 2 / 8.), (Rw <= 3).astype(float)
cnt = np.array([(lab == s).sum() for s in names], float)
sh = np.round(kr * n).astype(int); band = (kr / smpd <= 1. / a.lpmask * 1.5)   # shells kept for the regression: to 2/3 of LPMASK
shb = sh[band]; nshb = shb.max() + 1; wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[band]
EF = np.array([np.fft.rfftn(M * rd('%s/recvol_state%02d_even_unfil.mrc' % (a.rdir, s)))[band] for s in names]); OF = np.array([np.fft.rfftn(M * rd('%s/recvol_state%02d_odd_unfil.mrc' % (a.rdir, s)))[band] for s in names])
lpb = lpf[band]
def back(Fb):   # band coefficients -> low-passed real-space map inside the sphere
    x = np.zeros(kr.shape, complex); x[band] = Fb * lpb; return np.fft.irfftn(x, s=(n, n, n))[inm]
def shell_sum(z): return np.bincount(shb, weights=wb * z, minlength=nshb)
T = None
if a.gt and a.truth:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); gi = G[:, 0].astype(int) - 1; ok = (gi >= 0) & (gi < len(lab)); conf[gi[ok]] = G[ok, 1].astype(int)
    tc = {}; T = []
    for s in names:
        cs, cn = np.unique(conf[lab == s], return_counts=True); good = [(c, k) for c, k in zip(cs, cn) if c > 0]; acc = 0; tot = 0
        for c, k in good:
            if c not in tc:
                fs = [x for d in a.truth.split(',') for x in sorted(glob.glob('%s/conf%02d_*m160.mrc' % (d, c)))]; tc[c] = lp(M * rd(fs[0]))[inm] if fs else None
            if tc[c] is not None: acc = acc + k * tc[c]; tot += k
        T.append(acc / tot if tot else np.full(inm.sum(), np.nan))
    T = np.array(T); tc.clear()
def full(v): x = np.zeros((n, n, n)); x[inm] = v; return x
def smooth(v): return conv(full(v), gk)[inm]
rows = ['# node leaves particles | tail mask % of sphere, noise thr | reproducible PCs, PC1 share, PC1/PC2 support overlap | truth: recall of the true variance, true support % of sphere (90 % of the variance) | alternative masks: size, recall | child masks inside this one']
masks = {}
nodes = [r for r in sorted(mem) if len(mem[r]) >= a.minleaves]
for u in sorted(nodes, key=lambda r: -len(mem[r])):
    ms = mem[u]; w = cnt[ms] / cnt[ms].sum(); EuF, OuF = w @ EF[ms], w @ OF[ms]; den = shell_sum(np.real(EuF * OuF.conj())); okd = den > 1e-6 * den.max(); den = np.where(okd, den, 1.)
    def coef(z): return np.where(okd, z / den, 1.)[shb]   # per-shell regression coefficient on the node mean (1 where the mean has no signal)
    dE = np.array([back(EF[j] - coef(shell_sum(np.real(EF[j] * OuF.conj()))) * EuF) for j in ms]); dO = np.array([back(OF[j] - coef(shell_sum(np.real(OF[j] * EuF.conj()))) * OuF) for j in ms])
    v = smooth((w[:, None] * dE * dO).sum(0)); thr = np.quantile(-v, 1 - a.alpha); core = v > thr
    vp = np.maximum(v, 0); oo = np.argsort(-vp); cm = np.cumsum(vp[oo]) / max(vp.sum(), 1e-300); core80 = np.zeros(len(v), bool); core80[oo[:np.searchsorted(cm, 0.8) + 1]] = True; core10 = v > 0.1 * v.max()
    def soft(c): m_ = (conv(full(c.astype(float)), ball) > 0.5 / ball.sum()).astype(float); return np.clip(conv(m_, gk), 0, 1) * M
    vs_ = np.where(core, v, 0.); o2 = np.argsort(-vs_); c2 = np.cumsum(vs_[o2]) / max(vs_.sum(), 1e-300); sig80 = np.zeros(len(v), bool); sig80[o2[:np.searchsorted(c2, 0.8) + 1]] = True; sig80 &= core   # significant AND holding 80 % of the significant variance
    alt = {'mass80': soft(core80), 'rel10': soft(core10), 'sig80': soft(sig80)}
    mk = (conv(full(core.astype(float)), ball) > 0.5 / ball.sum()).astype(float); mk = np.clip(conv(mk, gk), 0, 1) * M; masks[u] = mk
    wr('%s/masks/node%03d.mrc' % (a.out, u), mk); wr('%s/masks/node%03d_sig80.mrc' % (a.out, u), alt['sig80']); fr = 100 * mk.sum() / M.sum()
    # cross-half PCA of the members' deviations (weights sqrt(w) so the eigen-decomposition is of the weighted covariance)
    sw = np.sqrt(w); Gm = (sw[:, None] * dE) @ (sw[:, None] * dO).T; Gm = (Gm + Gm.T) / 2; ev, U = np.linalg.eigh(Gm); ev, U = ev[::-1], U[:, ::-1]
    nrep = int((ev > max(-ev.min(), 0)).sum()); sh1 = ev[0] / max(ev[ev > 0].sum(), 1e-300)
    ov = np.nan
    if nrep >= 2:
        pcs = [smooth(((sw * U[:, k])[:, None] * dE).sum(0) * ((sw * U[:, k])[:, None] * dO).sum(0)) for k in range(2)]
        sup = [p > np.quantile(-p, 1 - a.alpha) for p in pcs]; ov = (sup[0] & sup[1]).sum() / max(min(sup[0].sum(), sup[1].sum()), 1)
    tr = ''
    if T is not None and np.isfinite(T[ms]).all():
        Tu = w @ T[ms]; tv = smooth((w[:, None] * (T[ms] - Tu) ** 2).sum(0)); tv = np.maximum(tv, 0); o = np.argsort(-tv); cs_ = np.cumsum(tv[o]) / max(tv.sum(), 1e-300)
        tsup = np.zeros(len(tv), bool); tsup[o[:np.searchsorted(cs_, 0.9) + 1]] = True; mkin = mk[inm]
        tr = '%.2f %.1f' % ((mkin * tv).sum() / max(tv.sum(), 1e-300), 100 * tsup.sum() / inm.sum()) + ' | ' + ' '.join('%s %.1f%% rec %.2f' % (k_, 100 * m_.sum() / M.sum(), (m_[inm] * tv).sum() / max(tv.sum(), 1e-300)) for k_, m_ in alt.items())
    inside = ''
    if u in kids:
        inside = ' '.join('%d:%.2f' % (c, (masks.get(c, np.zeros(1)) * mk).sum() / max(masks.get(c, np.zeros(1)).sum(), 1e-9)) for c in kids[u] if c in masks)
    rows.append('%d %d %d | %.1f %.3g | %d %.2f %s | %s | %s' % (u, len(ms), cnt[ms].sum(), fr, thr, nrep, sh1, '-' if np.isnan(ov) else '%.2f' % ov, tr or '-', inside or '-'))
# child-inside-parent: a second pass now that every mask exists
out = [rows[0]]
for l in rows[1:]:
    u = int(l.split()[0]); w_ = l.rsplit('|', 1)
    if u in kids: w_[1] = ' ' + (' '.join('%d:%.2f' % (c, (masks[c] * masks[u]).sum() / max(masks[c].sum(), 1e-9)) for c in kids[u] if c in masks) or '-')
    out.append('|'.join(w_))
open(a.out + '/node_masks.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
