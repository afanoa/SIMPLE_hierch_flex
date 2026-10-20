#!/usr/bin/env python3
"""Is a pair difference a spectral (radial filter) difference of the pooled map? Per shell: cross-tree FSC of the raw
difference, correlation of each tree's raw difference with its pooled map (|corr| ~ 1 = a pure per-shell scale/filter
difference), and the cross-tree FSC after removing a per-shell multiple of the pooled map. usage: pairdump_spectral.py NPZ"""
import sys, numpy as np
z = np.load(sys.argv[1]); c = z['rawA'].shape[0]; f = np.fft.fftfreq(c) * c / 160.; kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * 160).astype(int), 80).ravel()
env = z['env']; FA, FB = np.fft.fftn(z['rawA'] * env).ravel(), np.fft.fftn(z['rawB'] * env).ravel(); PA = np.fft.fftn(z['PA'] * env).ravel()
def shcorr(X, Y):
    num = np.bincount(sh, np.real(X * np.conj(Y)), 81); dx = np.bincount(sh, np.abs(X) ** 2, 81); dy = np.bincount(sh, np.abs(Y) ** 2, 81); return num / np.maximum(np.sqrt(dx * dy), 1e-30)
def deproj(X, P):   # remove the per-shell least-squares multiple of P
    a = np.bincount(sh, np.real(X * np.conj(P)), 81) / np.maximum(np.bincount(sh, np.abs(P) ** 2, 81), 1e-30); return X - a[sh] * P
print('pair statistic %.2f' % float(z['stat']))
f0 = shcorr(FA, FB); cA = shcorr(FA, PA); cB = shcorr(FB, PA); f1 = shcorr(deproj(FA, PA), deproj(FB, PA))
print('%5s %7s %9s %9s %9s %12s' % ('shell', 'res A', 'FSC raw', 'corr(DA,P)', 'corr(DB,P)', 'FSC deproj'))
for k in list(range(2, 16)) + list(range(16, 41, 4)): print('%5d %7.1f %9.3f %9.3f %9.3f %12.3f' % (k, 160 * 1.464 / k, f0[k], cA[k], cB[k], f1[k]))
