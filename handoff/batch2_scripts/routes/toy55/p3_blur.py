#!/usr/bin/env python3
"""Local nuisance regression before the focus mask (problem 3). At realistic SNR the members of a node differ reproducibly over the whole
molecule while their true differences sit in a fraction of it, so something other than conformation varies between states; a candidate
is LOCAL blur / scale (states aligned with different accuracy are blurred differently, more at the periphery), which the per-shell
(global) regression of node_masks2/3 cannot remove. After that regression every member's half-deviation d is regressed locally (Gaussian
window WIN A) on the node mean m and its Laplacian L (both from the node mean, 20 A): variants none | blur (L) | scale (m) | scale+blur
(m, L); the residuals' cross-half variance -> the sig80 mask (as node_masks3). A motion is a gradient term (kept by all variants), an
occupancy change a local scale term (removed by 'scale'). Reports mask size (% sphere / % molecule), with truth the recall of the true
variance and the share of the cross-half variance outside the true support, with --arm the closed / open arm coverage.
usage: p3_blur.py OUT --list members.txt [--labels L --gt G --truth D] [--arm C:CL] [--name N --mskdiam 220 --win 15]"""
import os, glob, argparse, warnings, numpy as np
warnings.filterwarnings('ignore')
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--list', required=True); ap.add_argument('--labels', default=''); ap.add_argument('--gt', default='')
ap.add_argument('--truth', default=''); ap.add_argument('--arm', default=''); ap.add_argument('--name', default='set'); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--win', type=float, default=15.); ap.add_argument('--alpha', type=float, default=0.001); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
L_ = [l.split() for l in open(a.list) if l.strip() and not l.startswith('#')]; ids = [int(l[0]) for l in L_]; cnt = np.array([float(l[3]) for l in L_]); w = cnt / cnt.sum()
n = int(np.fromfile(L_[0][1], np.int32, 1)[0]); smpd = float(np.fromfile(L_[0][1], np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = min(a.mskdiam / 2 / smpd, n / 2 - 7)
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2); kA = kr / smpd
sh = np.round(kr * n).astype(int); nsh = sh.max() + 1; wg = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)
def ssum(z): return np.bincount(sh.ravel(), (wg * z).ravel(), nsh)
lpg = lambda r_: np.exp(-(kA * r_) ** 2 * np.log(2)); lp20 = lpg(20.)
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
def conv(x, k): return np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(k / k.sum()), s=(n, n, n))
gk, ball = np.exp(-Rw ** 2 / 8.), (Rw <= 3).astype(float); sw_ = a.win / smpd / 2.355; gwin = np.exp(-Rw ** 2 / (2 * sw_ ** 2)); Gwin = np.fft.rfftn(gwin / gwin.sum())
wconv = lambda x: np.fft.irfftn(np.fft.rfftn(x) * Gwin, s=(n, n, n))
def crop(x, m):
    N = x.shape[0]
    if N == m: return x
    F = np.fft.fftshift(np.fft.fftn(x)); lo = N // 2 - m // 2; return np.real(np.fft.ifftn(np.fft.ifftshift(F[lo:lo + m, lo:lo + m, lo:lo + m]))) * (m / N) ** 3
EF = [np.fft.rfftn(M * rd(l[1])) for l in L_]; OF = [np.fft.rfftn(M * rd(l[2])) for l in L_]
Eu, Ou = sum(wi * x for wi, x in zip(w, EF)), sum(wi * x for wi, x in zip(w, OF)); den = ssum(np.real(Eu * Ou.conj())); okd = np.abs(den) > 1e-6 * np.abs(den).max(); den = np.where(okd, den, 1.)
co = lambda x, y: np.where(okd, ssum(np.real(x * y.conj())) / den, 1.)[sh]
dE = [np.fft.irfftn((x - co(x, Ou) * Eu) * lp20, s=(n, n, n)) for x in EF]; dO = [np.fft.irfftn((x - co(x, Eu) * Ou) * lp20, s=(n, n, n)) for x in OF]
Mu = (Eu + Ou) / 2; m = np.fft.irfftn(Mu * lp20, s=(n, n, n)); Lp = np.fft.irfftn(-4 * np.pi ** 2 * kr ** 2 * Mu * lp20, s=(n, n, n))
mean10 = np.fft.irfftn(Mu * lpg(10.), s=(n, n, n)); supp = (mean10 > 0.2 * mean10[inm].max()) & inm
mC = mO = None
if a.arm:
    cons, clo = [np.fft.irfftn(np.fft.rfftn(crop(rd(p), n)) * lpg(10.), s=(n, n, n)) for p in a.arm.split(':')]; supc = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
    def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supc] ** 2)), 1e-12)
    d = std(clo) - std(cons); thr = 0.5 * np.percentile(np.abs(d[supc | (clo > 0.2 * clo.max())]), 99); mC, mO = d > thr, d < -thr
tv = None
if a.labels and a.gt and a.truth:
    lab = np.loadtxt(a.labels, dtype=int); G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); gi = G[:, 0].astype(int) - 1; ok = (gi >= 0) & (gi < len(lab)); conf[gi[ok]] = G[ok, 1].astype(int)
    tc = {}; TT = []
    for s in ids:
        cs, cn = np.unique(conf[lab == s], return_counts=True); acc, tot = 0., 0
        for c, k in zip(cs, cn):
            if c <= 0: continue
            if c not in tc:
                fs = sorted(glob.glob('%s/conf%02d_*m160.mrc' % (a.truth, c))); tc[c] = np.fft.irfftn(np.fft.rfftn(M * crop(rd(fs[0]), n)) * lp20, s=(n, n, n)) if fs else None
            if tc[c] is not None: acc = acc + k * tc[c]; tot += k
        TT.append(acc / tot)
    tc.clear(); Tu = sum(wi * t for wi, t in zip(w, TT)); tv = np.maximum(conv(sum(wi * (t - Tu) ** 2 for wi, t in zip(w, TT)), gk), 0) * M
    o = np.argsort(-tv.ravel()); c_ = np.cumsum(tv.ravel()[o]) / tv.sum(); tsup = np.zeros(tv.size, bool); tsup[o[:np.searchsorted(c_, 0.9) + 1]] = True; tsup = tsup.reshape(tv.shape)
def sig80(v):
    thr = np.quantile(-v[inm], 1 - a.alpha); core = (v > thr) & inm; vs = np.where(core, v, 0.).ravel(); o = np.argsort(-vs); c = np.cumsum(vs[o]) / max(vs.sum(), 1e-300)
    k = np.zeros(vs.size, bool); k[o[:np.searchsorted(c, 0.8) + 1]] = True; k = k.reshape(v.shape) & core
    mk = (conv(k.astype(float), ball) > 0.5 / ball.sum()).astype(float); return np.clip(conv(mk, gk), 0, 1) * M
Smm, SmL, SLL = wconv(m * m), wconv(m * Lp), wconv(Lp * Lp)
def resid(d, kind):
    if kind == 'none': return d
    if kind == 'blur': b = wconv(d * Lp) / np.maximum(SLL, 1e-12 * SLL.max()); return d - b * Lp
    if kind == 'scale': s_ = wconv(d * m) / np.maximum(Smm, 1e-12 * Smm.max()); return d - s_ * m
    p, q = wconv(d * m), wconv(d * Lp); det = np.maximum(Smm * SLL - SmL ** 2, 1e-12 * (Smm * SLL).max())
    s_ = (SLL * p - SmL * q) / det; b = (Smm * q - SmL * p) / det; return d - s_ * m - b * Lp
out = ['# %s (%d members, %d particles; molecule %.1f %% of the sphere, window %.0f A): variant | mask %%sphere/%%molecule | truth: recall, true support, off-support share of the cross-half variance (data / truth) | arm closed/open, share on arm'
       % (a.name, len(ids), cnt.sum(), 100. * supp.sum() / inm.sum(), a.win)]
if tv is not None: out.append('   true support (90 %% of the true variance): %.1f %%s / %.1f %%m' % (100. * tsup.sum() / inm.sum(), 100. * (tsup & supp).sum() / supp.sum()))
for kind in ('none', 'blur', 'scale', 'scale+blur'):
    rE, rO = [resid(d, kind) for d in dE], [resid(d, kind) for d in dO]; v = conv(sum(wi * x * y for wi, x, y in zip(w, rE, rO)), gk); mk = sig80(v); b = mk > 0.5
    s = '%-10s | %5.1f%%s / %5.1f%%m' % (kind, 100. * b.sum() / inm.sum(), 100. * (b & supp).sum() / supp.sum())
    if tv is not None:
        vp = np.maximum(v, 0) * M; s += ' | recall %.2f, off-support share %.2f' % ((mk * tv).sum() / tv.sum(), vp[~tsup].sum() / max(vp.sum(), 1e-300))
    if mC is not None: s += ' | arm %.2f/%.2f, on-arm %.2f' % ((b & mC).sum() / max(mC.sum(), 1), (b & mO).sum() / max(mO.sum(), 1), (b & (mC | mO)).sum() / max(b.sum(), 1))
    out.append(s); print(s, flush=True)
open('%s/%s_blur.txt' % (a.out, a.name), 'w').write('\n'.join(out) + '\n'); print('\n'.join(out[:2]))
