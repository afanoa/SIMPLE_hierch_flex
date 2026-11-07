#!/usr/bin/env python3
"""usage: quality_screen2.py CASES   (quality_screen.py's case list; the lower-corr group L is read from QSCREEN lines on stdin)
Map-only quality signatures of D = M_L - M_H (15 A low-pass, the case's mask; no particle corr involved):
- dB: the blur difference in the stable region (mask voxels with abs(D) under 3 sd of its noise): per shell from 15 to 6 A the
  log amplitude ratio |M_L| / |M_H| of the soft-masked stable region, fitted against s^2; dB = -4 x slope (A^2, > 0 = L blurrier),
  its SE from the fit;
- dens: Spearman of abs(D) with the mean map's density over the mask (a blur or amplitude change tracks density everywhere; a
  conformational change sits where the parts move);
- core: the share of the top 1 % abs(D) voxels inside the densest 20 % of the mask;
- rad: the median radius of the top 1 % abs(D) voxels over the mask's radius (0 = centre, 1 = edge).
One QS2 line per case."""
import sys, re, numpy as np
from scipy.stats import spearmanr
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = h.view(np.float32); n = int(h[0])
    return np.frombuffer(raw[1024 + int(h[23]):], np.float32, count=n ** 3).astype(np.float64).reshape(n, n, n), float(hf[10]) / n
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; G = F[s:s + n, s:s + n, s:s + n]
    else: G = np.zeros((n, n, n), complex); s = n // 2 - c; G[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / m) ** 3
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def soft(m, n): return np.real(np.fft.ifftn(np.fft.fftn(m.astype(float)) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2)))
LB = {}
for l in sys.stdin:
    m = re.match(r'QSCREEN (\S+)\s+\S+\s+\| L = ([AB])', l)
    if m: LB[m.group(1)] = m.group(2) == 'B'
for l in open(sys.argv[1]):
    w = l.split()
    if len(w) < 11 or w[0].startswith('#') or w[0] not in LB: continue
    MA, MB = [read_mrc(f)[0] for f in w[2:5]], [read_mrc(f)[0] for f in w[5:8]]; smpd = read_mrc(w[2])[1]; n = MA[0].shape[0]; ML, MH = (MB, MA) if LB[w[0]] else (MA, MB)
    mk = fresample(read_mrc(w[8])[0], n) > 0.5; K = kr(n); sh = np.rint(K).astype(int); res = 15.
    lpw = np.clip((1.1 * smpd / res - K / n) / (0.2 * smpd / res), 0., 1.); lpw = 0.5 - 0.5 * np.cos(np.pi * lpw)
    FL, FH = [np.fft.fftn(x) for x in ML], [np.fft.fftn(x) for x in MH]
    Dr = np.real(np.fft.ifftn((FL[0] - FH[0]) * lpw)); Nr = np.real(np.fft.ifftn(((FL[1] - FL[2]) - (FH[1] - FH[2])) * lpw)) / 2; sd = Nr[mk].std()
    Mm = np.real(np.fft.ifftn((FL[0] + FH[0]) / 2 * lpw)); stab = mk & (np.abs(Dr) <= 3 * sd)
    if stab.sum() < 0.1 * mk.sum(): stab = mk & (np.abs(Dr) <= np.percentile(np.abs(Dr[mk]), 30))
    ss = soft(stab, n); AL, AH = np.abs(np.fft.fftn(ML[0] * ss)), np.abs(np.fft.fftn(MH[0] * ss)); s1, s2 = int(np.ceil(n * smpd / 15.)), int(n * smpd / 6.); x, y = [], []
    for k in range(s1, s2 + 1):
        q = sh == k; x.append((k / (n * smpd)) ** 2); y.append(np.log(np.sqrt((AL[q] ** 2).sum() / (AH[q] ** 2).sum())))
    x, y = np.array(x), np.array(y); A = np.column_stack([np.ones_like(x), x]); c, r, *_ = np.linalg.lstsq(A, y, rcond=None); res_ = y - A @ c
    se = np.sqrt(res_ @ res_ / max(len(y) - 2, 1) / ((x - x.mean()) ** 2).sum()); dB = -4 * c[1]; dBse = 4 * se
    ad = np.abs(Dr[mk]); rho = spearmanr(ad, Mm[mk])[0]; top = ad >= np.percentile(ad, 99); dense = Mm[mk] >= np.percentile(Mm[mk], 80)
    g = (np.arange(n) - n // 2) * smpd; Z, Y, X = np.meshgrid(g, g, g, indexing='ij'); rr = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)[mk]
    print('QS2 %-36s %-8s | dB %+.1f A^2 (SE %.1f, z %+.1f) | dens rho %+.3f | top-1%% in densest 20%%: %.3f | rad %.2f' % (w[0], w[1], dB, dBse, dB / max(dBse, 1e-9), rho,
          (top & dense).sum() / top.sum(), np.median(rr[top]) / rr.max())); sys.stdout.flush()
