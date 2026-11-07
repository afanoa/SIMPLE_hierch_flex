#!/usr/bin/env python3
"""check_notch.py: FSC of two half maps per shell after Fourier-cropping both to box B, their shell power, and the ring power /
mean CTF^2 of a cached plane set. usage: check_notch.py EVEN ODD B [CACHE.npy]"""
import sys, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo')
from mrcio import read_mrc, fcrop3
ev, od, B = sys.argv[1], sys.argv[2], int(sys.argv[3])
a, smpd = read_mrc(ev); b, _ = read_mrc(od)
print('maps', a.shape, 'smpd', smpd, '-> box', B)
if a.shape[0] != B: a = fcrop3(a, B); b = fcrop3(b, B)
A = np.fft.fftn(np.fft.ifftshift(a)); Bm = np.fft.fftn(np.fft.ifftshift(b))
f = np.fft.fftfreq(B) * B
r = np.rint(np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)).astype(int)
for s in range(1, B // 2 + 1):
    m = r == s
    xy = (A[m] * np.conj(Bm[m])).real.sum(); xx = (np.abs(A[m]) ** 2).sum(); yy = (np.abs(Bm[m]) ** 2).sum()
    print('shell %2d res %6.1f fsc %7.3f power(even) %10.3e' % (s, B * a.shape[0] * smpd / a.shape[0] / max(s, 1), xy / max(np.sqrt(xx * yy), 1e-30), xx / m.sum()))
if len(sys.argv) > 4:
    c = np.load(sys.argv[4], mmap_mode='r'); n, Bc, nh = c.shape
    kk, hh = np.meshgrid(np.r_[0:Bc // 2, -Bc // 2:0], np.arange(nh), indexing='ij'); rr = np.rint(np.sqrt(hh ** 2 + kk ** 2)).astype(int)
    pw = np.zeros(Bc // 2 + 1); cnt = np.zeros(Bc // 2 + 1)
    for i in range(0, n, 10000):
        x = np.asarray(c[i:i + 2000]); p = (x.real ** 2 + x.imag ** 2).mean(0)
        pw += np.bincount(rr.ravel(), p.ravel(), Bc // 2 + 1); cnt += np.bincount(rr.ravel(), minlength=Bc // 2 + 1)
    print('plane ring power (mean |X|^2 per sample, 12k particles):')
    for s in range(1, Bc // 2 + 1): print('ring %2d %10.4e' % (s, pw[s] / max(cnt[s], 1)))
