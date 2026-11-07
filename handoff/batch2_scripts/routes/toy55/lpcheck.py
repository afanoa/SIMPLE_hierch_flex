#!/usr/bin/env python3
"""usage: lpcheck.py SITE_DIR .. -- for each site's bin1..5.mrc vs lp10/bin1..5.mrc: the share of Fourier power beyond 10 A and at
12-20 A kept by the filtered copy, and the real-space correlation of the two (a 10 A low-pass keeps the low band and removes the high one)"""
import sys, numpy as np
def rd(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = h.view(np.float32); n = int(h[0])
    return np.frombuffer(raw[1024 + int(h[23]):], np.float32, count=n ** 3).astype(float).reshape(n, n, n), float(hf[10]) / n
for d in sys.argv[1:]:
    for k in range(1, 6):
        a, s = rd('%s/bin%d.mrc' % (d, k)); b, _ = rd('%s/lp10/bin%d.mrc' % (d, k)); n = a.shape[0]; f = np.fft.fftfreq(n) / s
        r = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2); A = np.abs(np.fft.fftn(a)) ** 2; B = np.abs(np.fft.fftn(b)) ** 2
        hi = r > 1 / 10.; lo = (r > 1 / 20.) & (r < 1 / 12.)
        print('%s bin%d: power kept beyond 10 A %.4f, at 12-20 A %.3f; corr %.3f' % (d.split('/')[-1], k, B[hi].sum() / A[hi].sum(), B[lo].sum() / A[lo].sum(), np.corrcoef(a.ravel(), b.ravel())[0, 1]))
