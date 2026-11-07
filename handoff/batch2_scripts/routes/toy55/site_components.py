#!/usr/bin/env python3
"""Reproducible variability components of a node from its many small groups (no segmentation, no thresholds on the maps). Every group's
even / odd map is low-passed to LP A and scaled to the weighted mean inside the support; the halves' weighted deviations
X_E[g] = sqrt(w_g) (E_g - E), X_O[g] = sqrt(w_g) (O_g - O) give the K x K cross-half matrix M = X_E X_O^T, symmetrised; its
eigenvectors u are the patterns of change shared by the two independent halves (eigenvalue = reproducible variance), and the negative
eigenvalues, which only noise produces, give the null: components above the largest |negative eigenvalue| are reproducible. Per
component: the half maps X_E^T u, X_O^T u (their correlation = reproducibility), the map (their mean; OUT/comp_NN.mrc), its
localisation (share of the support holding 50 / 90 % of its energy) and the group loadings (a minority = few groups carry it:
loading kurtosis, share of |u|^2 in the top 3 groups). Optional truth (--truth name:A.mrc:B.mrc): per component the share of its
energy inside each truth region (|A - B| > 30 % of its max) and the correlation of the component with A - B.
usage: site_components.py OUT --list LIST (name even odd particles) [--lp 15 --mskdiam 220 --truth n:A:B ...]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--list', required=True); ap.add_argument('--lp', type=float, default=15.)
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--truth', nargs='*', default=[]); ap.add_argument('--nshow', type=int, default=8)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
L = [l.split() for l in open(a.list) if l.strip() and not l.startswith('#')]; w = np.array([float(l[3]) for l in L]); w /= w.sum(); K = len(L)
n = int(np.fromfile(L[0][1], np.int32, 1)[0]); smpd = float(np.fromfile(L[0][1], np.float32, 13)[10]) / n
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) / smpd
lpf = np.exp(-(k * a.lp) ** 2 * np.log(2)); lpr = lambda x: np.fft.irfftn(np.fft.rfftn(x) * lpf, s=x.shape)
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); sph = R <= a.mskdiam / 2 / smpd
E = np.array([lpr(rd(l[1]))[sph] for l in L]); O = np.array([lpr(rd(l[2]))[sph] for l in L]); Em, Om = w @ E, w @ O
mean_ = 0.5 * (Em + Om); sup = mean_ > 0.1 * mean_.max()
for j in range(K):
    E[j] *= (E[j][sup] @ Em[sup]) / max(E[j][sup] @ E[j][sup], 1e-30); O[j] *= (O[j][sup] @ Om[sup]) / max(O[j][sup] @ O[j][sup], 1e-30)
Em, Om = w @ E, w @ O; XE = np.sqrt(w)[:, None] * (E - Em); XO = np.sqrt(w)[:, None] * (O - Om)
M = XE @ XO.T; S = 0.5 * (M + M.T); lam, U = np.linalg.eigh(S); o = np.argsort(-lam); lam, U = lam[o], U[:, o]
null = max(abs(lam[lam < 0]).max() if (lam < 0).any() else 0., 1e-30); rep = lam > null
print('%d groups: cross-half eigenvalues (reproducible variance, x 1e3) %s; null (largest |negative|) %.3g -> %d reproducible components'
      % (K, ' '.join('%.3g' % (1e3 * v) for v in lam[:12]), 1e3 * null, rep.sum()))
TR = []
for t in a.truth:
    nm, fa, fb = t.split(':'); A_, B_ = rd(fa), rd(fb)
    if A_.shape[0] != n:
        def rs(x):
            nb = x.shape[0]; F = np.fft.fftshift(np.fft.fftn(x)); out = np.zeros((n, n, n), complex); c = min(nb, n); lo, lo2 = nb // 2 - c // 2, n // 2 - c // 2
            out[lo2:lo2 + c, lo2:lo2 + c, lo2:lo2 + c] = F[lo:lo + c, lo:lo + c, lo:lo + c]; return np.real(np.fft.ifftn(np.fft.ifftshift(out))) * (n / nb) ** 3
        A_, B_ = rs(A_), rs(B_)
    d = (lpr(A_) - lpr(B_))[sph]; TR.append((nm, d, np.abs(d) > 0.3 * np.abs(d).max()))
for i in range(min(a.nshow, K)):
    u = U[:, i]; cE, cO = XE.T @ u, XO.T @ u; r = cE @ cO / max(np.linalg.norm(cE) * np.linalg.norm(cO), 1e-30); c = 0.5 * (cE + cO)
    e = np.sort(c ** 2)[::-1]; ce = np.cumsum(e) / e.sum(); n50, n90 = np.searchsorted(ce, 0.5) + 1, np.searchsorted(ce, 0.9) + 1
    uu = u ** 2; top3 = np.sort(uu)[::-1][:3].sum(); kurt = ((u - u.mean()) ** 4).mean() / max(((u - u.mean()) ** 2).mean() ** 2, 1e-30)
    full = np.zeros((n, n, n)); full[sph] = c; wr(os.path.join(a.out, 'comp_%02d.mrc' % (i + 1)), full, smpd)
    tr = '; '.join('%s: energy share %.2f, corr %+.2f' % (nm, (c[m] ** 2).sum() / (c ** 2).sum(), (c @ d) / max(np.linalg.norm(c) * np.linalg.norm(d), 1e-30)) for nm, d, m in TR)
    print('   comp %d%s: eigenvalue %.2f x null, half-map corr %.2f, 50 / 90 %% of its energy in %.1f / %.1f %% of the support; loadings: top 3 groups %.2f of |u|^2 (kurtosis %.1f), largest %s%s'
          % (i + 1, '' if rep[i] else ' (null)', lam[i] / null, r, 100 * n50 / sup.sum(), 100 * n90 / sup.sum(), top3, kurt,
             ' '.join('%s:%+.2f' % (L[j][0], u[j]) for j in np.argsort(-np.abs(u))[:4]), ('; ' + tr) if tr else ''))
np.savetxt(os.path.join(a.out, 'loadings.txt'), U[:, :min(a.nshow, K)], fmt='%+.4f', header='group loadings (rows = groups in list order), columns = components by eigenvalue')
