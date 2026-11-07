#!/usr/bin/env python3
"""Tree-weighted FINAL maps of a run's delivered units (merged leaves). Unit half-maps = particle-weighted averages of the member leaves'
unfiltered half-maps (as unit_maps.py); a tree over the units (state_merge.py --list, complete linkage, in OUT/tree). Per Fourier shell the
units' cross-half covariance C_ab = Re<E_a, O_b> is fitted with the tree (one deviation power per edge plus the root, NNLS; tree_wiener5.py)
or taken as measured (PSD projection: emp), in two regions (inside / outside the focus mask = where the units differ between halves); each
unit's map = the posterior mean C (C + N)^-1 x (x = half-mean, N its noise; halves with the half-map noise 2N, each half borrowing from the
same half of the others). Estimators: raw (= unit_maps.py), own (each unit's own per-shell Wiener), tree2, emp2. The tree's off-diagonal R^2
over HP-LPFIT A (vs a shuffled tree) picks tree2 when >= R2MIN, else emp2 (OUT/choice.txt). Full maps padded in Fourier space to BOX,
raised-cosine low-pass LP, unsharpened -> <MAPDIR>_<est>/[junk_]unit_NN.mrc for the DELIVER estimators (OUT/maps_<est> for the others),
the unfiltered padded halves in its halves/ subfolder.
usage: tw_units.py OUT --stage halves --leaves leaves_int.txt --units final_leaves.txt [--junk final_junk.txt]      (unit half-maps + list)
       tw_units.py OUT --stage maps --mapdir BASE [--deliver tree2,emp2 --box 160 --lp 10 --mskdiam 232 --hp 40 --lpfit 8 --r2min 0.7]"""
import os, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import nnls
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--stage', required=True, choices=['halves', 'maps']); ap.add_argument('--leaves', default='')
ap.add_argument('--units', default=''); ap.add_argument('--junk', default=''); ap.add_argument('--mapdir', default=''); ap.add_argument('--deliver', default='tree2,emp2')
ap.add_argument('--box', type=int, default=160); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--mskdiam', type=float, default=232.)
ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lpfit', type=float, default=8.); ap.add_argument('--r2min', type=float, default=0.7)
ap.add_argument('--lpmask', type=float, default=20.); ap.add_argument('--mthr', type=float, default=0.15); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd_):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd_; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
if a.stage == 'halves':
    L = {int(l.split()[0]): (l.split()[1], l.split()[2], float(l.split()[3])) for l in open(a.leaves) if l.strip() and not l.startswith('#')}
    U = [(int(l.split()[0]), [int(x) for x in l.split()[1].split(',')]) for l in open(a.units) if l.strip() and not l.startswith('#')]
    os.makedirs(a.out + '/unit_halves', exist_ok=True); rows = []
    for u, mem in U:
        c = np.array([L[x][2] for x in mem]); w = c / c.sum(); smpd = float(np.fromfile(L[mem[0]][0], np.float32, 13)[10]) / int(np.fromfile(L[mem[0]][0], np.int32, 1)[0])
        for h, j in (('even', 0), ('odd', 1)): wr('%s/unit_halves/unit_%02d_%s.mrc' % (a.out, u, h), sum(wi * rd(L[x][j]) for wi, x in zip(w, mem)), smpd)
        rows.append('%d %s/unit_halves/unit_%02d_even.mrc %s/unit_halves/unit_%02d_odd.mrc %d' % (u, a.out, u, a.out, u, c.sum()))
    open(a.out + '/units_list.txt', 'w').write('\n'.join(rows) + '\n'); print('%d unit half-map pairs (box %d, %.3f A/px) from %d leaves' % (len(U), rd(L[U[0][1][0]][0]).shape[0], smpd, len(L))); raise SystemExit
# ---- stage maps: units (list), their tree
LST = [l.split() for l in open(a.out + '/units_list.txt') if l.strip()]; files = {int(w[0]): (w[1], w[2]) for w in LST}; ucnt = {int(w[0]): float(w[3]) for w in LST}
junk = {int(x) for x in open(a.junk).read().split()} if a.junk and os.path.isfile(a.junk) else set()
tf = a.out + '/tree/merge_tree.txt'; hdr = open(tf).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; NL = len(names)
S = np.atleast_2d(np.loadtxt(tf, comments='#')); par = {}
for st in S: par[int(st[1])] = int(st[3]); par[int(st[2])] = int(st[3])
root = int(S[-1, 3])
def path(i):
    p = []
    while i != root: p.append(i); i = par[i]
    return p
cnt = np.array([ucnt[s] for s in names])
f0 = files[names[0]][0]; n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
E = [M * rd(files[s][0]) for s in names]; O = [M * rd(files[s][1]) for s in names]
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
sh = np.round(kr * n).astype(int); sel = (sh >= 1) & (sh < n // 2); order = np.argsort(sh[sel], kind='stable'); shs = sh[sel][order]; nsh = n // 2
wg = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[sel][order]; bounds = np.searchsorted(shs, np.arange(1, nsh + 1)); starts = bounds[:-1]
def ft(x): return np.fft.rfftn(x)[sel][order]
def ift(v):
    full = np.zeros(kr.shape, complex); tmp = np.zeros(sel.sum(), complex); tmp[order] = v; full[sel] = tmp; return np.fft.irfftn(full, s=(n, n, n))
Ef = np.array([ft(x) for x in E]); Of = np.array([ft(x) for x in O])
# ---- focus mask: where the units differ between independent halves (as tree_wiener5)
wl = cnt / cnt.sum(); lpf = np.exp(-(kr / smpd * a.lpmask) ** 2 * np.log(2)); el = [np.fft.irfftn(np.fft.rfftn(x) * lpf, s=(n, n, n)) for x in E]
ol = [np.fft.irfftn(np.fft.rfftn(x) * lpf, s=(n, n, n)) for x in O]; eb = sum(w * x for w, x in zip(wl, el)); ob = sum(w * x for w, x in zip(wl, ol))
v = sum(w * (x - eb) * (y - ob) for w, x, y in zip(wl, el, ol)); del el, ol
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
vs = conv(v, np.exp(-Rw ** 2 / 8.)); core = (vs > a.mthr * vs[inm].max()) & inm; Mv = (conv(core.astype(float), (Rw <= 3).astype(float)) > 0.5 / (Rw <= 3).sum()).astype(float)
Mv = np.clip(conv(Mv, np.exp(-Rw ** 2 / 8.)), 0, 1) * M; wr(a.out + '/focus_mask.mrc', Mv, smpd)
EA = np.array([ft(x * Mv) for x in E]); OA = np.array([ft(x * Mv) for x in O]); del E, O
# ---- covariance structures (tree_wiener5)
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
    yo = y[offd]; return Cf, 1 - ((yo - yf[offd]) ** 2).sum() / max(((yo - yo.mean()) ** 2).sum(), 1e-300)
pick = {int(round(n * smpd / r)): r for r in (20., 10., 6.)}
def wiener_all(Ef_, Of_, keep=False):   # full-map (x) and half-map (e, o) posterior means per estimator
    X_ = (Ef_ + Of_) / 2; out = {k: [np.zeros_like(X_) for _ in range(3)] for k in ('own', 'tree', 'emp')}; Wk = {}; fit = {k: [] for k in ('star', 'tree', 'shuf')}
    for s in range(1, nsh):
        lo, hi = bounds[s - 1], bounds[s]; w = wg[lo:hi]; e, o, x = Ef_[:, lo:hi], Of_[:, lo:hi], X_[:, lo:hi]; sw = w.sum()
        Ce = np.real((e * w) @ o.conj().T + (o * w) @ e.conj().T) / (2 * sw); Nn = np.maximum((np.abs(e - o) ** 2 * w).sum(1) / (4 * sw), 1e-300)
        d = np.maximum(np.diag(Ce), 0); out['own'][0][:, lo:hi] = (d / (d + Nn))[:, None] * x
        out['own'][1][:, lo:hi] = (d / (d + 2 * Nn))[:, None] * e; out['own'][2][:, lo:hi] = (d / (d + 2 * Nn))[:, None] * o
        ev, U = np.linalg.eigh(Ce); Cs = {'emp': (U * np.maximum(ev, 0)) @ U.T}
        for k, A in (('star', Ast), ('tree', Atr), ('shuf', Ash)):
            Cf, r2 = covfit(A, Ce)
            if k == 'tree': Cs[k] = Cf
            if a.lpfit <= n * smpd / s <= a.hp: fit[k].append(r2)
        for k, C in Cs.items():
            W = np.linalg.solve(C + np.diag(Nn), C).T; Wh = np.linalg.solve(C + np.diag(2 * Nn), C).T
            out[k][0][:, lo:hi] = W @ x; out[k][1][:, lo:hi] = Wh @ e; out[k][2][:, lo:hi] = Wh @ o
            if keep and s in pick: Wk[(k, s)] = W
    return out, Wk, {k: float(np.mean(v_)) for k, v_ in fit.items()}
est, Wk, fit = wiener_all(Ef, Of, keep=True); oA, _, _ = wiener_all(EA, OA); oB, _, _ = wiener_all(Ef - EA, Of - OA)
EST = {'raw': [(Ef + Of) / 2, Ef, Of], 'own': est['own'], 'tree2': [oA['tree'][j] + oB['tree'][j] for j in range(3)], 'emp2': [oA['emp'][j] + oB['emp'][j] for j in range(3)]}
del oA, oB, EA, OA
choice = 'tree2' if fit['tree'] >= a.r2min else 'emp2'
# ---- per-unit half-map FSC of every estimator (precision readout; borrowed halves share the weights, not the noise)
def fsc_res(eh, oh):
    num = np.add.reduceat(wg * np.real(eh * oh.conj()), starts); den = np.sqrt(np.add.reduceat(wg * np.abs(eh) ** 2, starts) * np.add.reduceat(wg * np.abs(oh) ** 2, starts))
    fsc = num / np.maximum(den, 1e-300); b = np.where(fsc < 0.143)[0]; res_sh = n * smpd / np.arange(1, nsh); return res_sh[b[0] - 1] if len(b) and b[0] > 0 else res_sh[-1]
real = [i for i in range(NL) if names[i] not in junk]
out = ['%d units (%d junk), box %d at %.3f A/px -> delivered at box %d, low-pass %.0f A, unsharpened; focus mask %.1f %% of the sphere' % (NL, len(junk), n, smpd, a.box, a.lp, 100 * Mv.sum() / M.sum()),
       'covariance fit, off-diagonal R^2 over %.0f-%.0f A: tree %.3f, shuffled tree %.3f, star %.3f -> deliver %s (rule: tree2 when the tree R^2 >= %.2f, else the unstructured emp2)' % (a.hp, a.lpfit, fit['tree'], fit['shuf'], fit['star'], choice, a.r2min),
       'half-map FSC 0.143 resolution, median / worst over the %d non-junk units (precision only): %s' % (len(real), ' | '.join('%s %.1f / %.1f A' % (k, np.median([fsc_res(EST[k][1][i], EST[k][2][i]) for i in real]), max(fsc_res(EST[k][1][i], EST[k][2][i]) for i in real)) for k in ('raw', 'own', 'tree2', 'emp2')))]
i0 = min(real, key=lambda i: cnt[i]); dist = np.array([len(set(P_tree[i0]) ^ set(P_tree[j])) for j in range(NL)])
for k in ('tree', 'emp'):
    for s_, r in sorted(pick.items()):
        if (k, s_) in Wk:
            Wr = Wk[(k, s_)][i0]; out.append('  smallest unit %d (%d ptcls) at %4.1f A, %s: own weight %.2f, from the others %.2f (by tree distance %s)' % (names[i0], cnt[i0], r, k, Wr[i0], Wr.sum() - Wr[i0], ' '.join('d%d %.2f' % (d_, Wr[dist == d_].sum()) for d_ in np.unique(dist) if d_ > 0)))
# ---- write: padded to BOX, full maps low-passed (raised cosine, as unit_maps.py), halves unfiltered in halves/
nb = a.box; sb = smpd * n / nb if nb > n else smpd
fb = np.fft.fftfreq(max(nb, n)); qx, qy, qz = np.meshgrid(fb, fb, np.fft.rfftfreq(max(nb, n)), indexing='ij'); qk = np.sqrt(qx ** 2 + qy ** 2 + qz ** 2) / sb
k1, k2 = 1. / (a.lp + 1.), 1. / (a.lp - 0.5); filt = np.where(qk <= k1, 1., np.where(qk >= k2, 0., 0.5 + 0.5 * np.cos(np.pi * (qk - k1) / (k2 - k1)))) if a.lp > 0 else np.ones_like(qk)
def pad(V):
    if nb <= n: return V
    F = np.fft.fftshift(np.fft.fftn(V)); P = np.zeros((nb,) * 3, complex); lo = nb // 2 - n // 2; P[lo:lo + n, lo:lo + n, lo:lo + n] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(P))) * (nb / n) ** 3
deliver = [k for k in a.deliver.split(',') if k]; dirs = {}
for k in ('raw', 'own', 'tree2', 'emp2'):
    d_ = '%s_%s' % (a.mapdir, k) if k in deliver else '%s/maps_%s' % (a.out, k); os.makedirs(d_ + '/halves', exist_ok=True); dirs[k] = d_
    for i, s in enumerate(names):
        nm = '%sunit_%02d' % ('junk_' if s in junk else '', s); V = pad(ift(EST[k][0][i]))
        wr('%s/%s.mrc' % (d_, nm), np.fft.irfftn(np.fft.rfftn(V) * filt, s=V.shape), sb)
        for j, h in ((1, 'even'), (2, 'odd')): wr('%s/halves/%s_%s.mrc' % (d_, nm, h), pad(ift(EST[k][j][i])), sb)
open(a.out + '/choice.txt', 'w').write('%s %s\n' % (choice, dirs[choice])); open(a.out + '/map_dirs.txt', 'w').write(''.join('%s %s\n' % (k, d_) for k, d_ in dirs.items()))
open(a.out + '/summary.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out)); print('maps: ' + ' '.join('%s=%s' % kv for kv in dirs.items()))
