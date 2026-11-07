#!/usr/bin/env python3
"""Gold-standard check of independently discovered states: state maps of half A and half B (and optionally a reference
run, e.g. f692e on all particles) as deviations from their own run's mean state, 12 A low-pass, inside the consensus
support; A<->B mutual best matches above the 99th percentile of the non-matching correlations = reproduced states; a
reference state is supported when it matches a state in BOTH halves above that threshold. Closed-arm occupancy from each
run's arm_states.txt. usage: match_halves.py DIR_A DIR_B [DIR_REF]  (DIR = a consolidated_* dir with arm_states.txt)"""
import sys, os, glob, numpy as np
def rd(f):
    h = np.fromfile(f, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(f, np.float32, offset=1024 + int(h[23]), count=int(nx * ny * nz)).reshape(nz, ny, nx).astype(np.float64)
def maps_of(d):
    arm = {}
    for l in open(d + '/arm_states.txt'):
        if l.startswith('#') or not l.strip(): continue
        p = l.split(); arm[p[0]] = float(p[1])
    fs = sorted(arm); return fs, np.array([arm[f] for f in fs])
dirs = sys.argv[1:4]; runs = [maps_of(d) for d in dirs]
n = rd(runs[0][0][0]).shape[0]; fr = np.fft.fftfreq(n, d=1.464); kx, ky, kz = np.meshgrid(fr, fr, fr, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (12. / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
def devs(fs):
    V = np.array([lp(rd(f)) for f in fs]); V = (V - V.mean(axis=(1, 2, 3), keepdims=True)) / V.std(axis=(1, 2, 3), keepdims=True); return V, V.mean(0)
D, M = [], []
for fs, _ in runs: V, m = devs(fs); D.append(V); M.append(m)
supp = sum(M) / len(M); supp = supp > np.percentile(supp, 80)   # consensus support: the densest 20 % of the box
def vec(V, m): X = (V - m)[:, supp]; X = X - X.mean(1, keepdims=True); return X / np.linalg.norm(X, axis=1, keepdims=True)
X = [vec(V, m) for V, m in zip(D, M)]
C = X[0] @ X[1].T; ba, bb = C.argmax(1), C.argmax(0); mutual = [(i, ba[i]) for i in range(len(ba)) if bb[ba[i]] == i]
mask = np.ones_like(C, bool)
for i, j in mutual: mask[i, j] = False
thr = float(np.percentile(C[mask], 99)) if mask.any() else 0.
rep = [(i, j, C[i, j]) for i, j in mutual if C[i, j] >= thr]
print('half A %d states, half B %d states; mutual best pairs %d; threshold (99th pct of the other pairs) %.3f -> REPRODUCED %d' % (C.shape[0], C.shape[1], len(mutual), thr, len(rep)))
armA, armB = runs[0][1], runs[1][1]
for i, j, c in sorted(rep, key=lambda t: -t[2]): print('   A %-40s B %-40s corr %.3f  closed-arm A %.2f B %.2f' % (os.path.basename(os.path.dirname(os.path.dirname(runs[0][0][i]))) + '/' + os.path.basename(runs[0][0][i]), os.path.basename(runs[1][0][j]), c, armA[i], armB[j]))
print('closed-arm >= 0.8: half A %d, half B %d; reproduced pairs with both >= 0.8: %d, both >= 0.6: %d' % ((armA >= 0.8).sum(), (armB >= 0.8).sum(), sum(1 for i, j, c in rep if armA[i] >= 0.8 and armB[j] >= 0.8), sum(1 for i, j, c in rep if armA[i] >= 0.6 and armB[j] >= 0.6)))
if len(X) == 3:
    CA, CB = X[2] @ X[0].T, X[2] @ X[1].T; ok = (CA.max(1) >= thr) & (CB.max(1) >= thr); armR = runs[2][1]
    print('reference %d states: supported in BOTH halves %d (in A only %d, in B only %d, neither %d)' % (len(ok), ok.sum(), ((CA.max(1) >= thr) & ~(CB.max(1) >= thr)).sum(), (~(CA.max(1) >= thr) & (CB.max(1) >= thr)).sum(), (~(CA.max(1) >= thr) & ~(CB.max(1) >= thr)).sum()))
    for lo in (0.8, 0.6):
        sel = armR >= lo; print('   reference states with closed-arm >= %.1f: %d, supported in both halves %d' % (lo, sel.sum(), (ok & sel).sum()))
