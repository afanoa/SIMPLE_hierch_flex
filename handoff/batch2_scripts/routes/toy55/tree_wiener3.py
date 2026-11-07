#!/usr/bin/env python3
"""Tree-weighted reconstruction, post hoc on a refine3D_states run (map space): a test of borrowing along the state tree.
Model: every node's map = its parent's map + a deviation; per resolution shell one deviation power per tree edge plus the root's
signal power, fitted (non-negative least squares) to the leaves' cross-half signal covariances C_ab = Re <E_a, O_b> (unbiased: the
halves' noise is independent). Each leaf's map is the posterior mean V = C (C + N)^-1 x per shell (x = the leaves' half-means,
N = their noise from |E - O|^2): every leaf borrows from every other, the weights falling with tree distance and, per shell, where
relatives disagree. Compared against truth (each leaf's composition-weighted truth map, junk excluded) with: raw half-means, each
leaf's own Wiener filter, a star (shrinkage to the consensus), the tree with shuffled leaves, the unstructured covariance, and the
tree in two regions (inside / outside the focus mask = where the leaves differ between halves: their cross-half variance map).
Scores: FSC 0.5 resolution vs truth, real-space correlation vs truth, correlation of each leaf's deviation from the consensus with
its true deviation (HP-LP band); the tree's fit to the measured covariances (off-diagonal R^2, 40-10 A) needs no truth.
tree_wiener3.py: --gt/--truth optional (real data: the truth-free covariance fit and borrowing readout only); --maps K1,K2 writes
those estimators' leaf maps (OUT/maps/<est>/state_NN.mrc; raw = the half-means through the same masking) for readouts such as arm_readout.
usage: tree_wiener3.py OUT --rdir RD --labels refined_states.txt --tree merge_tree.txt --gt gt_particles.txt --truth dir1,dir2
       [--mskdiam 220 --lpmask 20 --mthr 0.15 --hp 40 --lp 8 --seed 1]"""
import os, glob, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import nnls
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', required=True); ap.add_argument('--labels', required=True)
ap.add_argument('--tree', required=True); ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--maps', default='')
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--lpmask', type=float, default=20.); ap.add_argument('--mthr', type=float, default=0.15)
ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lp', type=float, default=8.); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
lab = np.loadtxt(a.labels, dtype=int); conf = np.full(len(lab), -1); HT = bool(a.gt and a.truth)
if HT: G = np.loadtxt(a.gt, comments='#'); gi = G[:, 0].astype(int) - 1; ok = (gi >= 0) & (gi < len(lab)); conf[gi[ok]] = G[ok, 1].astype(int)
hdr = open(a.tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; NL = len(names)
S = np.atleast_2d(np.loadtxt(a.tree, comments='#')); par = {}
for st in S: par[int(st[1])] = int(st[3]); par[int(st[2])] = int(st[3])
root = int(S[-1, 3])
def path(i):   # the edges (child node ids) from leaf node i up to the root
    p = []
    while i != root: p.append(i); i = par[i]
    return p
# ---- maps: sphere-masked half-maps, truth per leaf
f0 = '%s/recvol_state%02d_even_unfil.mrc' % (a.rdir, names[0]); n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
tdirs = a.truth.split(',') if a.truth else []; tcache = {}
def tmap(c):
    if c not in tcache:
        fs = [x for d in tdirs for x in sorted(glob.glob('%s/conf%02d_*m160.mrc' % (d, c)))]; tcache[c] = rd(fs[0]) if fs else None
        if fs and tcache[c].shape != (n, n, n): raise SystemExit('truth map %s is not box %d' % (fs[0], n))
    return tcache[c]
E, O, T, cnt, jf = [], [], [], [], []
for s in names:
    m = lab == s; cnt.append(m.sum()); E.append(M * rd('%s/recvol_state%02d_even_unfil.mrc' % (a.rdir, s))); O.append(M * rd('%s/recvol_state%02d_odd_unfil.mrc' % (a.rdir, s)))
    if not HT: jf.append(0.); continue
    cs, cn = np.unique(conf[m], return_counts=True); good = [(c, k) for c, k in zip(cs, cn) if c > 0 and tmap(c) is not None]; ng = sum(k for _, k in good)
    jf.append(1 - ng / max(m.sum(), 1)); T.append(M * sum(k * tmap(c) for c, k in good) / ng if good else np.zeros((n, n, n)))
cnt = np.array(cnt, float); jf = np.array(jf); tcache.clear()
# ---- Fourier layout: rfft half-space coefficients of shells 1 .. n/2-1, sorted by shell
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
sh = np.round(kr * n).astype(int); sel = (sh >= 1) & (sh < n // 2); order = np.argsort(sh[sel], kind='stable'); shs = sh[sel][order]; nsh = n // 2
wg = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[sel][order]; bounds = np.searchsorted(shs, np.arange(1, nsh + 1)); starts = bounds[:-1]
res_sh = n * smpd / np.arange(1, nsh)   # resolution of shells 1 .. nsh-1
def ft(x): return np.fft.rfftn(x)[sel][order]
def ift(v):
    full = np.zeros(kr.shape, complex); tmp = np.zeros(sel.sum(), complex); tmp[order] = v; full[sel] = tmp; return np.fft.irfftn(full, s=(n, n, n))
Ef = np.array([ft(x) for x in E]); Of = np.array([ft(x) for x in O]); Tf = np.array([ft(x) for x in T]) if HT else None; X = (Ef + Of) / 2; del T
# ---- focus mask: where the leaves differ between independent halves (cross-half variance of the low-passed leaves about their mean)
wl = cnt / cnt.sum(); lpf = np.exp(-(kr / smpd * a.lpmask) ** 2 * np.log(2)); el = [np.fft.irfftn(np.fft.rfftn(x) * lpf, s=(n, n, n)) for x in E]
ol = [np.fft.irfftn(np.fft.rfftn(x) * lpf, s=(n, n, n)) for x in O]; eb = sum(w * x for w, x in zip(wl, el)); ob = sum(w * x for w, x in zip(wl, ol))
v = sum(w * (x - eb) * (y - ob) for w, x, y in zip(wl, el, ol)); del el, ol
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
vs = conv(v, np.exp(-Rw ** 2 / 8.)); core = (vs > a.mthr * vs[inm].max()) & inm; Mv = (conv(core.astype(float), (Rw <= 3).astype(float)) > 0.5 / (Rw <= 3).sum()).astype(float)
Mv = np.clip(conv(Mv, np.exp(-Rw ** 2 / 8.)), 0, 1) * M; wr(a.out + '/focus_mask.mrc', Mv); wr(a.out + '/xhalf_variance.mrc', vs * M)
EA = np.array([ft(x * Mv) for x in E]); OA = np.array([ft(x * Mv) for x in O]); del E, O
# ---- covariance structures
pairs = [(p, q) for p in range(NL) for q in range(p, NL)]; offd = np.array([p != q for p, q in pairs])
def design(paths):
    ids = sorted(set(u for p in paths for u in p)); col = {u: j + 1 for j, u in enumerate(ids)}; A = np.zeros((len(pairs), 1 + len(ids)))
    for r, (p, q) in enumerate(pairs):
        A[r, 0] = 1
        for u in set(paths[p]) & set(paths[q]): A[r, col[u]] = 1
    return A
P_tree = [path(i) for i in range(NL)]; perm = np.random.default_rng(a.seed).permutation(NL)
Atr = design(P_tree); Ash = design([P_tree[perm[i]] for i in range(NL)]); Ast = design([[i] for i in range(NL)])
def covfit(A, Ce):
    y = np.array([Ce[p, q] for p, q in pairs]); th, _ = nnls(A, y); yf = A @ th; Cf = np.zeros((NL, NL))
    for r, (p, q) in enumerate(pairs): Cf[p, q] = Cf[q, p] = yf[r]
    yo = y[offd]; r2 = 1 - ((yo - yf[offd]) ** 2).sum() / max(((yo - yo.mean()) ** 2).sum(), 1e-300)
    return Cf, r2
pick = {int(round(n * smpd / r)): r for r in (20., 10., 6.)}
def wiener_all(Ef_, Of_, keep=False):
    X_ = (Ef_ + Of_) / 2; out = {k: np.zeros_like(X_) for k in ('own', 'star', 'tree', 'shuf', 'emp')}; Wk = {}; fit = {k: [] for k in ('star', 'tree', 'shuf')}
    for s in range(1, nsh):
        lo, hi = bounds[s - 1], bounds[s]; w = wg[lo:hi]; e, o, x = Ef_[:, lo:hi], Of_[:, lo:hi], X_[:, lo:hi]; sw = w.sum()
        Ce = np.real((e * w) @ o.conj().T + (o * w) @ e.conj().T) / (2 * sw); Nn = np.maximum((np.abs(e - o) ** 2 * w).sum(1) / (4 * sw), 1e-300)
        d = np.maximum(np.diag(Ce), 0); out['own'][:, lo:hi] = (d / (d + Nn))[:, None] * x
        ev, U = np.linalg.eigh(Ce); Cs = {'emp': (U * np.maximum(ev, 0)) @ U.T}
        for k, A in (('star', Ast), ('tree', Atr), ('shuf', Ash)):
            Cs[k], r2 = covfit(A, Ce)
            if a.lp <= n * smpd / s <= a.hp: fit[k].append(r2)
        for k, C in Cs.items():
            Wm = np.linalg.solve(C + np.diag(Nn), C).T; out[k][:, lo:hi] = Wm @ x
            if keep and s in pick: Wk[(k, s)] = Wm
    return out, Wk, {k: float(np.mean(v_)) for k, v_ in fit.items()}
est, Wk, fit = wiener_all(Ef, Of, keep=True); oA, _, _ = wiener_all(EA, OA); oB, _, _ = wiener_all(Ef - EA, Of - OA)
est['raw'] = X; est['tree2'] = oA['tree'] + oB['tree']; est['star2'] = oA['star'] + oB['star']; est['emp2'] = oA['emp'] + oB['emp']; del oA, oB, EA, OA
if a.maps:
    for k in a.maps.split(','):
        os.makedirs('%s/maps/%s' % (a.out, k), exist_ok=True)
        for i in range(NL): wr('%s/maps/%s/state_%02d.mrc' % (a.out, k, names[i]), ift(est[k][i]))
if not HT:
    i0 = int(np.argmin(cnt)); dist = np.array([len(set(P_tree[i0]) ^ set(P_tree[j])) for j in range(NL)])
    out = ['%d leaves, box %d, %.3f A/px; focus mask %.1f %% of the sphere' % (NL, n, smpd, 100 * Mv.sum() / M.sum()),
           'covariance fit, off-diagonal R^2 over %.0f-%.0f A (truth-free test of the tree): tree %.3f, shuffled tree %.3f, star %.3f' % (a.hp, a.lp, fit['tree'], fit['shuf'], fit['star'])]
    for k in ('tree', 'star'):
        for s_, r in sorted(pick.items()):
            if (k, s_) in Wk:
                Wr = Wk[(k, s_)][i0]; pp = Wr / cnt; pp /= pp[i0]
                out.append('  smallest leaf (%d, %d ptcls) at %4.1f A, %s: weight sum by tree distance %s | per particle, relative to its own: %s' % (names[i0], cnt[i0], r, k, ' '.join('d%d %.2f' % (d, Wr[dist == d].sum()) for d in np.unique(dist)), ' '.join('d%d %.2f' % (d, pp[dist == d].mean()) for d in np.unique(dist))))
    open(a.out + '/summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out)); raise SystemExit
# ---- scores against truth
tr = [ift(t) for t in Tf]; bl, bh = np.searchsorted(-res_sh, -a.hp), np.searchsorted(-res_sh, -a.lp, side='right'); blo, bhi = bounds[bl], bounds[bh]
Xb = (cnt[:, None] * X).sum(0) / cnt.sum(); Tb = (cnt[:, None] * Tf).sum(0) / cnt.sum()
lo20 = bounds[int(round(n * smpd / 20.))]
def scores(V):
    num = np.add.reduceat(wg * np.real(V * Tf.conj()), starts, axis=1); den = np.sqrt(np.add.reduceat(wg * np.abs(V) ** 2, starts, axis=1) * np.add.reduceat(wg * np.abs(Tf) ** 2, starts, axis=1))
    fsc = num / np.maximum(den, 1e-300); rs = []
    for i in range(NL):
        b = np.where(fsc[i] < 0.5)[0]; rs.append(res_sh[b[0] - 1] if len(b) and b[0] > 0 else (res_sh[0] if len(b) else res_sh[-1]))
    rcc = [np.corrcoef(ift(V[i])[inm], tr[i][inm])[0, 1] for i in range(NL)]; dv = []
    for i in range(NL):
        d1, d2 = (V[i] - Xb)[blo:bhi], (Tf[i] - Tb)[blo:bhi]; w = wg[blo:bhi]; e2 = (w * np.abs(d2) ** 2).sum()
        dv.append(np.nan if e2 < 1e-6 * (w * np.abs(Tf[i, blo:bhi]) ** 2).sum() else (w * np.real(d1 * d2.conj())).sum() / np.sqrt((w * np.abs(d1) ** 2).sum() * e2))
    return np.array(rs), np.array(rcc), np.array(dv)
fr = [(w_ * np.real(X[i, :lo20] * Tf[i, :lo20].conj())).sum() / np.sqrt((w_ * np.abs(X[i, :lo20]) ** 2).sum() * (w_ * np.abs(Tf[i, :lo20]) ** 2).sum()) for w_ in [wg[:lo20]] for i in range(NL)]
keys = ['raw', 'own', 'star', 'tree', 'tree2', 'star2', 'shuf', 'emp', 'emp2']; sc = {k: scores(est[k]) for k in keys}
okl = jf < 0.5; small = [i for i in np.argsort(cnt) if okl[i]][:5]
out = ['%d leaves, box %d, %.3f A/px; frame check (corr to truth to 20 A) median %.3f min %.3f; focus mask %.1f %% of the sphere; leaves with junk >= 50 %%: %d' % (NL, n, smpd, np.nanmedian(np.array(fr)[okl]), np.nanmin(np.array(fr)[okl]), 100 * Mv.sum() / M.sum(), (~okl).sum()),
       'covariance fit, off-diagonal R^2 over %.0f-%.0f A (truth-free test of the tree): tree %.3f, shuffled tree %.3f, star %.3f' % (a.hp, a.lp, fit['tree'], fit['shuf'], fit['star']),
       'estimator: FSC0.5 res vs truth (all / 5 smallest), real-space corr vs truth, deviation corr %.0f-%.0f A' % (a.hp, a.lp)]
nm = lambda x: np.nanmean(x) if np.isfinite(x).any() else np.nan
for k in keys:
    r_, c_, d_ = sc[k]; out.append('  %-5s res %5.1f / %5.1f A   rcc %.3f / %.3f   dev %.3f / %.3f' % (k, r_[okl].mean(), r_[small].mean(), c_[okl].mean(), c_[small].mean(), nm(d_[okl]), nm(d_[small])))
i0 = small[0] if small else int(np.argmin(cnt)); dist = np.array([len(set(P_tree[i0]) ^ set(P_tree[j])) for j in range(NL)])
for k in ('tree', 'star'):
    for s_, r in sorted(pick.items()):
        if (k, s_) in Wk:
            Wr = Wk[(k, s_)][i0]; pp = Wr / cnt; pp /= pp[i0]   # weight per particle of every leaf, relative to the leaf's own particles
            out.append('  smallest leaf (%d, %d ptcls) at %4.1f A, %s: weight sum by tree distance %s | per particle, relative to its own: %s' % (names[i0], cnt[i0], r, k, ' '.join('d%d %.2f' % (d, Wr[dist == d].sum()) for d in np.unique(dist)), ' '.join('d%d %.2f' % (d, pp[dist == d].mean()) for d in np.unique(dist))))
open(a.out + '/summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
with open(a.out + '/per_leaf.txt', 'w') as fo:
    fo.write('# leaf n junk | FSC0.5 res vs truth: %s | dev corr: %s\n' % (' '.join(keys), ' '.join(keys)))
    for i in np.argsort(cnt): fo.write('%d %d %.2f | %s | %s\n' % (names[i], cnt[i], jf[i], ' '.join('%.1f' % sc[k][0][i] for k in keys), ' '.join('%.3f' % sc[k][2][i] for k in keys)))
