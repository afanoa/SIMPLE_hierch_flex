#!/usr/bin/env python3
"""Can focus masks localize on real data (problem 3)? For every node of a state tree (or one node of --list members) the member-variance
mask of node_masks3.py (members regressed per Fourier shell on the node mean, cross-half variance, significant voxels holding 80 % of it)
and three variants that could localize where it stays broad:
  band  - the same variance from band-passed deviations (40-20, 20-12, 12-8 A) instead of the 20 A low-pass;
  comp  - one mask per reproducible cross-half component of the members' deviations (eigenvectors of the members' cross-half Gram,
          eigenvalue above the magnitude of the most negative one): v_k = c_k^E c_k^O, same significance rule;
  reg   - every member rigidly registered to the node mean first (6 parameters fitted on 15 A maps at box 64, applied to both halves),
          testing whether state-specific misregistration spreads the differences over the whole molecule.
Masks are reported as % of the molecule (node mean at 10 A above 0.2 of its maximum) and, with --arm CONSENSUS:CLOSED (real RNase P), the
share of the closed-arm / open-arm regions they cover (closed conformer minus consensus, as nm3_job.sbatch) and the share of the mask on
the arm; with --gt/--truth (synthetic) the recall of the members' true variance.
usage: p3_masks.py OUT --rdir RD --labels L --tree T | --list FILE  [--name N --minleaves 6 --mskdiam 232 --arm C:CL --gt G --truth D --reg 1 (largest nodes registered)]"""
import os, glob, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
from scipy.ndimage import affine_transform
from scipy.optimize import minimize
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', default=''); ap.add_argument('--labels', default=''); ap.add_argument('--tree', default='')
ap.add_argument('--list', default=''); ap.add_argument('--name', default='set'); ap.add_argument('--minleaves', type=int, default=6); ap.add_argument('--mskdiam', type=float, default=232.)
ap.add_argument('--arm', default=''); ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--alpha', type=float, default=0.001); ap.add_argument('--reg', type=int, default=1)   # registration for the REG largest nodes
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
if a.list:
    L = [l.split() for l in open(a.list) if l.strip() and not l.startswith('#')]; names = [int(l[0]) for l in L]; files = {int(l[0]): (l[1], l[2]) for l in L}
    cnt = np.array([float(l[3]) for l in L]); NL = len(names); mem = {0: list(range(NL))}; kids = {}; nodes = [0]; lab = None
else:
    lab = np.loadtxt(a.labels, dtype=int); hdr = open(a.tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; NL = len(names)
    files = {s: ('%s/recvol_state%02d_even_unfil.mrc' % (a.rdir, s), '%s/recvol_state%02d_odd_unfil.mrc' % (a.rdir, s)) for s in names}
    cnt = np.array([(lab == s).sum() for s in names], float); S = np.atleast_2d(np.loadtxt(a.tree, comments='#')); mem = {i: [i] for i in range(NL)}; kids = {}
    for st in S: p, q, r = int(st[1]), int(st[2]), int(st[3]); mem[r] = mem[p] + mem[q]; kids[r] = (p, q)
    nodes = sorted([u for u in mem if len(mem[u]) >= a.minleaves], key=lambda u: -len(mem[u]))
f0 = files[names[0]][0]; n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = min(a.mskdiam / 2 / smpd, n / 2 - 7)
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2); kA = kr / smpd   # cycles / A
sh = np.round(kr * n).astype(int); nsh = sh.max() + 1; wg = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)
def ssum(z): return np.bincount(sh.ravel(), (wg * z).ravel(), nsh)
lpg = lambda r_: np.exp(-(kA * r_) ** 2 * np.log(2))            # Gaussian low-pass, half power at r_ A
def bpf(lo, hi): return np.clip(lpg(hi) - lpg(lo), 0, None)       # band-pass between lo (coarse) and hi (fine) A
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
gk, ball = np.exp(-Rw ** 2 / 8.), (Rw <= 3).astype(float)
def back(F, filt): return np.fft.irfftn(F * filt, s=(n, n, n))
def crop(x, m):   # Fourier crop of a box-N map to box m (numpy unnormalised FFT: scale (m/N)^3)
    N = x.shape[0]
    if N == m: return x
    F = np.fft.fftshift(np.fft.fftn(x)); lo = N // 2 - m // 2; return np.real(np.fft.ifftn(np.fft.ifftshift(F[lo:lo + m, lo:lo + m, lo:lo + m]))) * (m / N) ** 3
EF = {s: np.fft.rfftn(M * rd(files[s][0])) for s in names}; OF = {s: np.fft.rfftn(M * rd(files[s][1])) for s in names}
# arm regions (real RNase P) at this box
mC = mO = None
if a.arm:
    lp10 = lpg(10.); cons, clo = [np.fft.irfftn(np.fft.rfftn(crop(rd(p), n)) * lp10, s=(n, n, n)) for p in a.arm.split(':')]
    supc = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
    def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supc] ** 2)), 1e-12)
    d = std(clo) - std(cons); thr = 0.5 * np.percentile(np.abs(d[supc | (clo > 0.2 * clo.max())]), 99); mC, mO = d > thr, d < -thr
# truth per member (synthetic)
T = None
if a.gt and a.truth and lab is not None:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); gi = G[:, 0].astype(int) - 1; ok = (gi >= 0) & (gi < len(lab)); conf[gi[ok]] = G[ok, 1].astype(int); tc = {}; T = {}
    for s in names:
        cs, cn = np.unique(conf[lab == s], return_counts=True); acc, tot = 0., 0
        for c, k in zip(cs, cn):
            if c <= 0: continue
            if c not in tc:
                fs = sorted(glob.glob('%s/conf%02d_*m160.mrc' % (a.truth, c))); tc[c] = np.fft.rfftn(M * crop(rd(fs[0]), n)) if fs else None
            if tc[c] is not None: acc = acc + k * tc[c]; tot += k
        T[s] = acc / tot if tot else None
    tc.clear()
def sig80(v):   # significant (above the noise quantile of the negative side) AND holding 80 % of the significant variance -> soft mask
    thr = np.quantile(-v[inm], 1 - a.alpha); core = (v > thr) & inm; vs = np.where(core, v, 0.).ravel(); o = np.argsort(-vs); c = np.cumsum(vs[o]) / max(vs.sum(), 1e-300)
    k = np.zeros(vs.size, bool); k[o[:np.searchsorted(c, 0.8) + 1]] = True; k = k.reshape(v.shape) & core
    m = (conv(k.astype(float), ball) > 0.5 / ball.sum()).astype(float); return np.clip(conv(m, gk), 0, 1) * M
def describe(mk, supp, tv=None):
    b = mk > 0.5; s = '%5.1f%%' % (100. * (b & supp).sum() / max(supp.sum(), 1))
    if mC is not None: s += ' arm %.2f/%.2f on-arm %.2f' % ((b & mC).sum() / max(mC.sum(), 1), (b & mO).sum() / max(mO.sum(), 1), (b & (mC | mO)).sum() / max(b.sum(), 1))
    if tv is not None: s += ' rec %.2f' % ((mk * tv).sum() / max(tv.sum(), 1e-300))
    return s
# rigid registration of a member to the node mean (fit at box 64 on 15 A maps)
nd_ = min(64, n); lp15 = lpg(15.)
def rotm(p):
    t = np.linalg.norm(p)
    if t < 1e-12: return np.eye(3)
    k = p / t; K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]); return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
def xform(v, p, order):   # rotate about the box centre by rotation vector p[:3] (rad), then shift p[3:] (voxels of v)
    c = np.full(3, v.shape[0] // 2, float); Ri = rotm(p[:3]).T; return affine_transform(v, Ri, offset=c - Ri @ (c + p[3:]), order=order, mode='constant')
def register(F, ref_ds):
    mov = crop(np.fft.irfftn(F * lp15, s=(n, n, n)), nd_); rm = ref_ds - ref_ds.mean()
    cost = lambda p: -np.sum((xform(mov, p, 1) - mov.mean()) * rm) / np.sqrt(max(np.sum((xform(mov, p, 1) - mov.mean()) ** 2) * np.sum(rm ** 2), 1e-300))
    res = minimize(cost, np.zeros(6), method='Powell', options={'xtol': 1e-3, 'ftol': 1e-6, 'maxfev': 600}); p = res.x.copy(); p[3:] *= n / nd_   # shift in full-box voxels
    return p
rows = ['# %s: node leaves particles | lp20 member mask (%% molecule, closed/open arm covered, share on arm%s) | bands 40-20 / 20-12 / 12-8 | components: n_rep; top 3 (share: mask) | registered: |shift| A median/max, rot deg median/max -> lp20 mask'
        % (a.name, ', truth recall' if T is not None else '')]
for iu, u in enumerate(nodes):
    ms = [names[i] for i in mem[u]]; w = cnt[mem[u]] / cnt[mem[u]].sum()
    def deviations(E, O):   # members minus their per-shell cross-half regression on the node mean (Fourier)
        Eu, Ou = sum(wi * E[s] for wi, s in zip(w, ms)), sum(wi * O[s] for wi, s in zip(w, ms)); den = ssum(np.real(Eu * Ou.conj())); okd = np.abs(den) > 1e-6 * np.abs(den).max(); den = np.where(okd, den, 1.)
        co = lambda x, y: np.where(okd, ssum(np.real(x * y.conj())) / den, 1.)[sh]
        return {s: E[s] - co(E[s], Ou) * Eu for s in ms}, {s: O[s] - co(O[s], Eu) * Ou for s in ms}, Eu, Ou
    dE, dO, Eu, Ou = deviations(EF, OF); mean10 = back((Eu + Ou) / 2, lpg(10.)); supp = (mean10 > 0.2 * mean10[inm].max()) & inm
    tv = None
    if T is not None and all(T[s] is not None for s in ms):
        Tu = sum(wi * T[s] for wi, s in zip(w, ms)); tv = np.maximum(conv(sum(wi * back(T[s] - Tu, lpg(20.)) ** 2 for wi, s in zip(w, ms)), gk), 0) * M
    var = lambda filt, dE_, dO_: conv(sum(wi * back(dE_[s], filt) * back(dO_[s], filt) for wi, s in zip(w, ms)), gk)
    lp20 = lpg(20.); m0 = sig80(var(lp20, dE, dO)); out = [describe(m0, supp, tv)]
    out.append(' / '.join(describe(sig80(var(bpf(lo, hi), dE, dO)), supp, None) for lo, hi in ((40., 20.), (20., 12.), (12., 8.))))
    # components of the members' cross-half Gram (20 A low-pass)
    Wl = wg * lp20 ** 2; sw = np.sqrt(w); Gm = np.array([[sw[i] * sw[j] * np.sum(Wl * np.real(dE[si] * dO[sj].conj())) for j, sj in enumerate(ms)] for i, si in enumerate(ms)]); Gm = (Gm + Gm.T) / 2
    ev, U = np.linalg.eigh(Gm); ev, U = ev[::-1], U[:, ::-1]; nrep = int((ev > max(-ev.min(), 0)).sum()); pos = ev[ev > 0].sum(); cm = []
    for k in range(min(3, max(nrep, 1))):
        cE = sum(sw[i] * U[i, k] * dE[s] for i, s in enumerate(ms)); cO = sum(sw[i] * U[i, k] * dO[s] for i, s in enumerate(ms)); mk = sig80(conv(back(cE, lp20) * back(cO, lp20), gk)); cm.append(mk)
        out.append('c%d %.2f: %s' % (k + 1, ev[k] / max(pos, 1e-300), describe(mk, supp, None)))
    comp = ' | '.join(out[2:]); out = out[:2]
    if len(cm) >= 2:
        b1, b2 = cm[0] > 0.5, cm[1] > 0.5; comp += ' | c1/c2 overlap %.2f' % ((b1 & b2).sum() / max(min(b1.sum(), b2.sum()), 1))
    if len(cm): comp += ' | union %s' % describe(np.maximum.reduce(cm), supp, tv)
    reg = '-'
    if iu < a.reg:
        ref_ds = crop(back((Eu + Ou) / 2, lp15), nd_); P_ = {s: register((EF[s] + OF[s]) / 2, ref_ds) for s in ms}
        sh_A = np.array([np.linalg.norm(P_[s][3:]) * smpd for s in ms]); rot_d = np.array([np.rad2deg(np.linalg.norm(P_[s][:3])) for s in ms])
        ER = {s: np.fft.rfftn(M * xform(np.fft.irfftn(EF[s], s=(n, n, n)), P_[s], 3)) for s in ms}; OR = {s: np.fft.rfftn(M * xform(np.fft.irfftn(OF[s], s=(n, n, n)), P_[s], 3)) for s in ms}
        dEr, dOr, _, _ = deviations(ER, OR); reg = '%.2f/%.2f A, %.2f/%.2f deg -> %s' % (np.median(sh_A), sh_A.max(), np.median(rot_d), rot_d.max(), describe(sig80(var(lp20, dEr, dOr)), supp, tv))
        del ER, OR, dEr, dOr
    rows.append('%s node %d %d %d | %s | %s | n_rep %d; %s | %s' % (a.name, u, len(ms), cnt[mem[u]].sum(), out[0], out[1], nrep, comp, reg)); print(rows[-1], flush=True)
open('%s/%s_masks.txt' % (a.out, a.name), 'w').write('\n'.join(rows) + '\n')
