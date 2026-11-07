#!/usr/bin/env python3
"""usage: motion_score.py SERIES_DIR MASK MSKDIAM_A
For a motion series in SERIES_DIR (bin1..5.mrc, halves/binK_even.mrc / binK_odd.mrc):
- per bin, the FSC 0.5 / 0.143 resolution of its halves inside a soft sphere of MSKDIAM_A (cosine edge 6 px);
- the correlation of neighbouring bins, and of the two end bins, inside the tree's root MASK (Fourier-resampled to the map box),
  unfiltered and low-passed at 15 / 10 A."""
import sys, numpy as np
D, MK, MD = sys.argv[1], sys.argv[2], float(sys.argv[3])
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = np.frombuffer(raw[:1024], np.float32)
    nx, ny, nz, mode, nsym = int(h[0]), int(h[1]), int(h[2]), int(h[3]), int(h[23])
    v = np.frombuffer(raw[1024 + nsym:], {0: np.int8, 1: np.int16, 2: np.float32, 6: np.uint16}[mode], count=nx * ny * nz).astype(np.float64)
    return v.reshape(nz, ny, nx), float(hf[10]) / nx
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; G = F[s:s + n, s:s + n, s:s + n]
    else: G = np.zeros((n, n, n), complex); s = n // 2 - c; G[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / m) ** 3
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def lowpass(v, smpd, res):
    n = v.shape[0]; r = kr(n) / n; fc = smpd / res; w = np.clip((1.1 * fc - r) / (0.2 * fc), 0., 1.); w = 0.5 - 0.5 * np.cos(np.pi * w)
    return np.real(np.fft.ifftn(np.fft.fftn(v) * w))
def fsc(a, b, smpd):
    n = a.shape[0]; x = np.arange(n) - n // 2; r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); R = MD / 2 / smpd
    m = np.clip((R + 6 - r) / 6, 0, 1); m = 0.5 - 0.5 * np.cos(np.pi * m); A = np.fft.fftn(a * m); B = np.fft.fftn(b * m); sh = np.rint(kr(n)).astype(int)
    num = np.bincount(sh.ravel(), (A * B.conj()).real.ravel()); da = np.bincount(sh.ravel(), (np.abs(A) ** 2).ravel()); db = np.bincount(sh.ravel(), (np.abs(B) ** 2).ravel())
    c = num / np.sqrt(np.maximum(da * db, 1e-30)); out = []
    for t in (0.5, 0.143):
        k = next((i for i in range(1, n // 2) if c[i] < t), n // 2); out.append(n * smpd / k)
    return out
maps = []; print('bin  FSC 0.5 / 0.143 (A)')
for k in range(1, 6):
    v, smpd = read_mrc('%s/bin%d.mrc' % (D, k)); e, _ = read_mrc('%s/halves/bin%d_even.mrc' % (D, k)); o, _ = read_mrc('%s/halves/bin%d_odd.mrc' % (D, k))
    r5, r1 = fsc(e, o, smpd); maps.append(v); print('bin%d %6.2f / %6.2f   (box %d, %.3f A/px)' % (k, r5, r1, v.shape[0], smpd))
mask, _ = read_mrc(MK); mk = fresample(mask, maps[0].shape[0]) > 0.5
def cc(a, b): a = a[mk] - a[mk].mean(); b = b[mk] - b[mk].mean(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))
for lab, res in (('unfiltered', None), ('15 A', 15.), ('10 A', 10.)):
    V = [m if res is None else lowpass(m, smpd, res) for m in maps]
    print('correlation %-10s neighbours %s | ends (bin1 vs bin5) %.3f' % (lab, ' '.join('%d-%d %.3f' % (k + 1, k + 2, cc(V[k], V[k + 1])) for k in range(4)), cc(V[0], V[4])))
