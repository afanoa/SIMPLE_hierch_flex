#!/usr/bin/env python3
"""Focus-mask sources against truth at one tree node (user, 2026-10-21: PC-derived masks may be noisy; real-space variation of
volumes instead?). Truth T(x): weighted per-voxel variance of the low-passed truth maps over the node members' conformers (junk left
out). Estimated variance maps, same low-pass: 'states' = cross-half real-space variance of the node's flex STATE half-maps after
per-shell regression on their weighted mean (the focus_mask.py statistic, i.e. what the current masks threshold); 'kids' = the same
for the cheap view-balanced children's half-maps; 'pcs' = sum_k lambda_k PC_k(x)^2 of the flex principal-component maps (unit norm).
Masks: the node's current mask.mrc, focus_mask.py on the children's half-maps, and the PC variance at the same 80 % mass rule (+3 px
dilation). Scores: corr(v, T) inside the sphere; per mask size (% of the sphere), recall (share of the true variance inside) and
enrichment (recall / size); the truth's own 80 % mass support for reference.
usage: nodemask_check.py --node NODEDIR --kids CHEAPKIDS_DIR --rows ROWS --gt GT --tmaps DIR [--lp 20] [--mskdiam 220] [--tag NAME]"""
import argparse, glob, os, subprocess, sys, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--node', required=True); ap.add_argument('--kids', required=True); ap.add_argument('--rows', required=True)
ap.add_argument('--gt', required=True); ap.add_argument('--tmaps', required=True); ap.add_argument('--lp', type=float, default=20.); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--tag', default=''); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
mem = [l.split() for l in open(os.path.join(a.node, 'members.txt')) if l.strip()]
n = int(np.fromfile(mem[0][1], np.int32, 1)[0]); smpd = float(np.fromfile(mem[0][1], np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.)); inm = M > 0.5
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
band = kr / smpd <= 1.5 / a.lp; shb = np.round(kr * n).astype(int)[band]; nb = shb.max() + 1; wb = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.)[band]
lpb = np.exp(-(kr / smpd * a.lp) ** 2 * np.log(2))[band]
def back(Fb): x = np.zeros(kr.shape, complex); x[band] = Fb * lpb; return np.fft.irfftn(x, s=(n, n, n))[inm]
def lowp(x): return back(np.fft.rfftn(M * x)[band])
def ssum(z): return np.bincount(shb, weights=wb * z, minlength=nb)
def upsample(x):  # flex PC maps live at the flex box (same physical extent): centred Fourier padding to the node box
    m = x.shape[0]
    if m == n: return x
    F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x))); P = np.zeros((n, n, n), complex); o = n // 2 - m // 2; P[o:o + m, o:o + m, o:o + m] = F
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(P))))
def crossvar(L):  # focus_mask.py statistic: L = [(even, odd, count)]
    w = np.array([c for _, _, c in L], float); w /= w.sum(); EF = np.array([np.fft.rfftn(M * rd(e))[band] for e, _, _ in L]); OF = np.array([np.fft.rfftn(M * rd(o))[band] for _, o, _ in L])
    EuF, OuF = w @ EF, w @ OF; den = ssum(np.real(EuF * OuF.conj())); okd = den > 1e-6 * den.max(); den = np.where(okd, den, 1.)
    coef = lambda z: np.where(okd, z / den, 1.)[shb]; v = np.zeros(inm.sum())
    for j in range(len(L)): v += w[j] * back(EF[j] - coef(ssum(np.real(EF[j] * OuF.conj()))) * EuF) * back(OF[j] - coef(ssum(np.real(OF[j] * EuF.conj()))) * OuF)
    return v
# truth: node members' conformer mix
rows = np.loadtxt(a.rows).astype(int); G = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(len(rows), int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
cs, cc = np.unique(conf[(rows > 0) & (conf > 0)], return_counts=True); S1 = 0.; S2 = 0.; W = 0.
for c_, k_ in zip(cs, cc):
    fl = glob.glob(os.path.join(a.tmaps, ('conf%02d_*_m160.mrc' % c_) if c_ < 100 else ('conf%d_*_m160.mrc' % c_)))
    if not fl: continue
    x = lowp(rd(fl[0])); S1 = S1 + k_ * x; S2 = S2 + k_ * x * x; W += k_
T = np.maximum(S2 / W - (S1 / W) ** 2, 0.)
V = {'states (current source)': crossvar([(l[1], l[2], float(l[3])) for l in mem])}
lab = np.loadtxt(os.path.join(a.kids, 'children.txt')).astype(int); RR = sorted(glob.glob(os.path.join(a.kids, 'rec', '*_reconstruct3D')))[-1]; K = lab.max()
kl = [(os.path.join(RR, 'recvol_state%02d_even_unfil.mrc' % k), os.path.join(RR, 'recvol_state%02d_odd_unfil.mrc' % k), float((lab == k).sum())) for k in range(1, K + 1)]
V['cheap children (PC1 bins)'] = crossvar(kl)
FD = sorted(glob.glob(os.path.join(a.node, '*_flex_pca')))[-1]; ev = np.loadtxt(os.path.join(FD, 'flex_pca_eigenvalues.txt'), comments='#', ndmin=2)
vp = 0.
for k_, lam in ev:
    pf = os.path.join(FD, 'flex_pca_pc%03d.mrc' % int(k_))
    if os.path.exists(pf): x = lowp(upsample(rd(pf))); vp = vp + lam * x * x / max((x * x).sum(), 1e-30)
V['PC maps (sum lambda PC^2)'] = vp
def mass_mask(v, mass=0.8):   # top voxels holding MASS of the positive variance, +3 px dilation (focus_mask.py rule without the significance test)
    vs = np.maximum(v, 0); o = np.argsort(-vs); c = np.cumsum(vs[o]) / max(vs.sum(), 1e-300); keep = np.zeros(len(v), bool); keep[o[:np.searchsorted(c, mass) + 1]] = True
    full = np.zeros((n, n, n)); full[inm] = keep; gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
    ball = (Rw <= 3).astype(float); d = np.fft.irfftn(np.fft.rfftn(full) * np.fft.rfftn(ball), s=(n, n, n)) > 0.5; return d[inm]
tot = T.sum(); tm = mass_mask(T)
print('%snode %s: %d members, %d conformers; LP %.0f A; the true variance\'s own 80 %% mass support covers %.1f %% of the sphere'
      % ((a.tag + ' ') if a.tag else '', os.path.basename(a.node.rstrip('/')), (rows > 0).sum(), len(cs), a.lp, 100 * tm.mean()))
for nm, v in V.items(): print('  corr(estimated variance, true variance) %-28s %.3f' % (nm, np.corrcoef(v, T)[0, 1]))
ksk = os.path.join(a.kids, 'kids_list.txt'); open(ksk, 'w').write(''.join('%d %s %s %d\n' % (i + 1, e, o, c) for i, (e, o, c) in enumerate(kl)))
km = os.path.join(a.kids, 'kids_mask.mrc'); subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'focus_mask.py'), km, '--list', ksk, '--mskdiam', str(a.mskdiam)], capture_output=True)
masks = {'current mask.mrc (flex states)': rd(os.path.join(a.node, 'mask.mrc'))[inm] > 0.5, 'focus_mask.py on cheap children': rd(km)[inm] > 0.5 if os.path.exists(km) else None,
         'PC variance, 80 % mass': mass_mask(V['PC maps (sum lambda PC^2)'])}
for nm, m in masks.items():
    if m is None: print('  mask %-34s missing' % nm); continue
    rec = T[m].sum() / tot; print('  mask %-34s size %5.1f %% of the sphere, holds %5.1f %% of the true variance, enrichment %.2f' % (nm, 100 * m.mean(), 100 * rec, rec / max(m.mean(), 1e-9)))
