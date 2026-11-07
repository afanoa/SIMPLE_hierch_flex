#!/usr/bin/env python3
"""Run-to-run comparison of two klem runs on the same particles, orientation-free: cell volumes matched by the correlation
of their deviations from the run's mean cell, the main axis of A traced through B's grid, per-particle agreement through
the matching (soft cross-run cell-deviation correlation vs a particle-permuted null) and junk agreement.
usage: klem_match.py KDIR_A KDIR_B [--n1 7 --n2 7]"""
import sys, glob, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('a'); ap.add_argument('b'); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); a = ap.parse_args()
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz, mode = h[0], h[1], h[2], h[3]; assert mode == 2
    return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
def load(d):
    fs = sorted(glob.glob(d + '/cells/cell_*.mrc')); return [read_mrc(f) for f in fs]
def common(A, B):
    """runs at different bands have coarse cubes of different sizes: both cropped in Fourier space to the smaller cube's band"""
    n = min(A[0].shape[0], B[0].shape[0]); F = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
    def crop(v):
        if v.shape[0] == n: return v.ravel()
        c = v.shape[0] // 2; h = n // 2; X = F(v)[c - h:c - h + n, c - h:c - h + n, c - h:c - h + n]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X)))).ravel() * (n / v.shape[0]) ** 3
    return np.array([crop(v) for v in A]), np.array([crop(v) for v in B])
def corr(X, Y):
    X = X - X.mean(1, keepdims=True); Y = Y - Y.mean(1, keepdims=True); return (X @ Y.T) / np.sqrt((X * X).sum(1)[:, None] * (Y * Y).sum(1)[None, :])
nc = a.n1 * a.n2; A, B = common(load(a.a)[:nc], load(a.b)[:nc]); dA, dB = A - A.mean(0), B - B.mean(0); CAB = corr(dA, dB); iu = np.triu_indices(nc, 1)
print('mean-cell correlation between runs %.3f; deviation rms A %.3g B %.3g; within-run cell-deviation corr median A %.2f B %.2f' % (corr(A.mean(0)[None], B.mean(0)[None])[0, 0], dA.std(), dB.std(), np.median(corr(dA, dA)[iu]), np.median(corr(dB, dB)[iu])))
best = CAB.argmax(1); print('between-run best-match deviation corr: median %.2f (min %.2f, max %.2f); distinct B cells matched %d/%d' % (np.median(CAB.max(1)), CAB.max(1).min(), CAB.max(1).max(), len(set(best)), nc))
i1, i2 = np.arange(nc) // a.n2, np.arange(nc) % a.n2; bi1, bi2 = i1[best], i2[best]
for row in range(a.n1): print('  A row z1=%d -> B cells (z1,z2): %s' % (row, ' '.join('(%d,%d)' % (bi1[k], bi2[k]) for k in np.where(i1 == row)[0])))
RA, RB = np.load(a.a + '/klem_resp.npz'), np.load(a.b + '/klem_resp.npz'); pa, pb = RA['pinds'], RB['pinds']; common, ia, ib = np.intersect1d(pa, pb, return_indices=True)
rA, rB = RA['r'][ia].astype(float), RB['r'][ib].astype(float); onA, onB = rA[:, :nc].sum(1) > 0.5, rB[:, :nc].sum(1) > 0.5; m = onA & onB
amA, amB = rA[m, :nc].argmax(1), rB[m, :nc].argmax(1)
print('%d common particles, %d on the grid in both: argmax(B) == bestmatch(argmax(A)) %.3f (chance %.3f); within one cell of the match %.3f; junk agreement %.3f' % (
    len(common), m.sum(), np.mean(amB == best[amA]), 1 / nc, np.mean((np.abs(bi1[amA] - i1[amB]) <= 1) & (np.abs(bi2[amA] - i2[amB]) <= 1)), np.mean(onA == onB)))
soft = np.sum(rA[m, :nc] @ CAB * rB[m, :nc], 1); null = np.sum(rA[m, :nc][np.random.default_rng(0).permutation(m.sum())] @ CAB * rB[m, :nc], 1)
print('soft per-particle cross-run cell-deviation correlation: mean %.3f (particles permuted: %.3f)' % (soft.mean(), null.mean()))
