#!/usr/bin/env python3
"""usage: deposited_loc.py MASK MAP ..   (the deposited conformers, docked; e.g. maps160/docked_all/conf*_m160.mrc)
The yardstick for a 'minor' split: split_effect.py's local change (D = Mi - Mj at shells to 15 A in the soft MASK, less its per-shell
projection on the pair's mean map, over the mean map's power; noise-free maps, so plain powers) for every pair of deposited
conformers, and for each conformer's nearest neighbour (highest 15 A CC in the mask). Prints the distributions."""
import sys, numpy as np
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
fs = sys.argv[2:]; M0, smpd = read_mrc(fs[0]); n = M0.shape[0]; mk = (fresample(read_mrc(sys.argv[1])[0], n) > 0.5).astype(float)
sm = np.real(np.fft.ifftn(np.fft.fftn(mk) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2))); sh = np.rint(kr(n)).astype(int).ravel(); smax = int(n * smpd / 15.); sel = (sh >= 1) & (sh <= smax); shs = sh[sel]
Q = [np.fft.fftn(read_mrc(f)[0] * sm).ravel()[sel] for f in fs]; k = len(Q); L = np.full((k, k), np.nan); CC = np.full((k, k), np.nan)
for i in range(k):
    for j in range(i + 1, k):
        d, pm = Q[i] - Q[j], (Q[i] + Q[j]) / 2; num = np.bincount(shs, (d.conj() * pm).real, smax + 1); den = np.bincount(shs, np.abs(pm) ** 2, smax + 1)
        c = (num / np.maximum(den, 1e-30))[shs]; dl = d - c * pm; L[i, j] = L[j, i] = (np.abs(dl) ** 2).sum() / (np.abs(pm) ** 2).sum()
        CC[i, j] = CC[j, i] = (Q[i].conj() * Q[j]).real.sum() / np.sqrt((np.abs(Q[i]) ** 2).sum() * (np.abs(Q[j]) ** 2).sum())
nn = [L[i, np.nanargmax(CC[i])] for i in range(k)]; allp = L[np.triu_indices(k, 1)]
q = lambda x: ' '.join('%.4f' % v for v in np.percentile(x, [0, 25, 50, 75, 100]))
print('%d deposited conformers, %d pairs; local change at 15 A (min q1 median q3 max):' % (k, len(allp)))
print('  nearest neighbour (highest CC): %s' % q(nn)); print('  all pairs:                      %s' % q(allp))
