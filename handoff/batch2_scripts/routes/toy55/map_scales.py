#!/usr/bin/env python3
"""Per state map: least-squares scale against the particle-weighted mean of the given states, and the share of the difference that
scale + a Gaussian blur explain (band 40-10 A, inside a spherical mask). usage: map_scales.py RDIR LABELS STATE,STATE,... [--smpd 1.464]"""
import sys, numpy as np
rd_, labf, sts = sys.argv[1:4]; smpd = 1.464; lab = np.loadtxt(labf, dtype=int); S = [int(x) for x in sts.split(',')]
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(float)
M = {s: rd('%s/recvol_state%02d.mrc' % (rd_, s)) for s in S}; n = M[S[0]].shape[0]; w = np.array([(lab == s).sum() for s in S], float); w /= w.sum()
f = np.fft.fftfreq(n); k = np.sqrt(sum(g ** 2 for g in np.meshgrid(f, f, f, indexing='ij'))); res = np.where(k > 0, smpd / np.maximum(k, 1e-9), 1e9); band = (res <= 40) & (res >= 10)
x = np.arange(n) - n / 2; msk = (np.sqrt(sum(g ** 2 for g in np.meshgrid(x, x, x, indexing='ij'))) < 75).astype(float)
F = {s: np.fft.fftn(M[s] * msk) for s in S}; Fm = sum(wi * F[s] for wi, s in zip(w, S))
for s in S:
    a = np.real(np.vdot(Fm[band], F[s][band])) / np.real(np.vdot(Fm[band], Fm[band])); d0 = np.sum(np.abs(F[s][band] - Fm[band]) ** 2); d1 = np.sum(np.abs(F[s][band] - a * Fm[band]) ** 2)
    best = d1; bb = 0.
    for B in np.linspace(-200, 200, 81):   # scale x exp(-B k^2 / 4) blur, k in 1/A
        g = np.exp(-B * (k[band] / smpd) ** 2 / 4.); A = np.real(np.vdot(Fm[band] * g, F[s][band])) / np.real(np.vdot(Fm[band] * g, Fm[band] * g)); d = np.sum(np.abs(F[s][band] - A * g * Fm[band]) ** 2)
        if d < best: best, bb = d, B
    print('state %2d (%5d particles): scale %.3f; difference to the mean explained by scale %.0f %%, by scale + blur (B %+.0f) %.0f %%' % (s, (lab == s).sum(), a, 100 * (1 - d1 / d0), bb, 100 * (1 - best / d0)))
