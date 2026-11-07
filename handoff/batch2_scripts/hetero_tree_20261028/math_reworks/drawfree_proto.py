#!/usr/bin/env python3
"""The draw-free pair statistic on a node's saved stage-B half maps (full N per state and half), as the Fortran defines it:
per shell the cross-half product of the difference maps (signal), each state's half-difference power over 1/nA + 1/nB (per-
particle noise), the expected cross-half FSC at a common n, 3-shell smoothed, peak over the hold-out band. The maps are
Fourier-cropped to the split test's box (2 px per 10 A) and masked with the node mask; the per-half state counts come from
the split log. Prints the per-shell curve at the test's n and at several n to show the statistic's n dependence.
usage: drawfree_proto.py NODE_DIR [--lp 10] [--band klo khi] [--n N]"""
import sys, os, re, math, glob, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mask_disc_proto import read_mrc

def to_box(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2; h = n // 2
    if m > n: G = F[c - h:c - h + n, c - h:c - h + n, c - h:c - h + n]
    else:
        G = np.zeros((n, n, n), complex); cc = n // 2; hh = m // 2
        G[cc - hh:cc - hh + m, cc - hh:cc - hh + m, cc - hh:cc - hh + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))).astype(np.float32) * (n / m)**3

ap = argparse.ArgumentParser(); ap.add_argument('node'); ap.add_argument('--lp', type=float, default=10.)
ap.add_argument('--n', type=int, default=0); ap.add_argument('--band', type=int, nargs=2, default=None)
a = ap.parse_args(); d = a.node.rstrip('/')
evens = sorted(glob.glob(os.path.join(d, 'stageB', 'recvol_state??_even_unfil.mrc')))
states = [int(re.search(r'state(\d\d)', e).group(1)) for e in evens]
log = open(os.path.join(d, 'split', 'hetero_step.log'), errors='replace').read()
m = re.search(r'A: (\d+) node particles \(labels ([0-9 ]+)\); B: (\d+) \(labels ([0-9 ]+)\)', log)
nA = [int(x) for x in m.group(2).split()]; nB = [int(x) for x in m.group(4).split()]
mb = re.search(r're-test band (\d+)\.\.(\d+)', log) or re.search(r'test band (\d+)\.\.(\d+)', log)
band = tuple(a.band) if a.band else (int(mb.group(1)), int(mb.group(2)))
mneq = re.search(r'\[cross-half EM groups\] equal N (\d+)', log) or re.search(r'equal N after the EM (\d+)', log) or re.search(r'equal N (\d+) per group', log)
n_common = a.n or int(mneq.group(1))
told = re.search(r're-test: weakest pair ([0-9.]+) at z ([0-9.-]+)', log)
# maps at the test's box
E0, smpd = read_mrc(evens[0]); box0 = E0.shape[0]
boxt = min(box0, max(26, 2 * math.ceil(box0 * smpd / a.lp)))
mask, _ = read_mrc(os.path.join(d, 'mask.mrc')); mask = to_box(mask, boxt) if mask.shape[0] != boxt else mask
mask = np.clip(mask, 0, 1)
def load(fn):
    v, _ = read_mrc(fn); v = to_box(v, boxt) if v.shape[0] != boxt else v
    return np.fft.fftn(v * mask)
A = {s: load(e) for s, e in zip(states, evens)}; B = {s: load(e.replace('_even_', '_odd_')) for s, e in zip(states, evens)}
f = np.fft.fftfreq(boxt) * boxt; fx, fy, fz = np.meshgrid(f, f, f, indexing='ij'); sh = np.rint(np.sqrt(fx**2 + fy**2 + fz**2)).astype(int)
nyq = boxt // 2; sh = np.minimum(sh, nyq)
def shell_sum(x):
    return np.bincount(sh.ravel(), weights=x.ravel(), minlength=nyq + 1)[:nyq + 1]
s2 = {}
for s in states:
    i = states.index(s)
    s2[s] = shell_sum(np.abs(A[s] - B[s])**2) / (1. / nA[i] + 1. / nB[i])
def T_df(p, q, n):
    c = shell_sum(np.real((A[p] - A[q]) * np.conj(B[p] - B[q])))
    nn = (s2[p] + s2[q]) / n
    fcurve = c / np.maximum(np.maximum(c, 0) + nn, 1e-30)
    sm = np.array([fcurve[k:k + 3].mean() for k in range(band[0], max(band[1] - 2, band[0]) + 1)])
    return sm.max(), fcurve
def fsc_actual(p, q):
    da = A[p] - A[q]; db = B[p] - B[q]
    num = shell_sum(np.real(da * np.conj(db))); den = np.sqrt(shell_sum(np.abs(da)**2) * shell_sum(np.abs(db)**2))
    fcurve = num / np.maximum(den, 1e-30)
    sm = np.array([fcurve[k:k + 3].mean() for k in range(band[0], max(band[1] - 2, band[0]) + 1)])
    return sm.max(), fcurve
print('node %s: box %d -> %d, %d states, per-half counts A %s B %s, band %d..%d, common n %d, old one-draw T %s (z %s)' % (
    d, box0, boxt, len(states), nA, nB, band[0], band[1], n_common, told.group(1) if told else '-', told.group(2) if told else '-'))
for i, p in enumerate(states):
    for q in states[i + 1:]:
        t, fc = T_df(p, q, n_common); ta, fa = fsc_actual(p, q)
        print('  pair %d-%d: draw-free T at n=%d: %.3f | full-N actual cross-half FSC peak: %.3f | T at n = %s: %s' % (
            p, q, n_common, t, ta, ' '.join(str(x) for x in [n_common // 4, n_common // 2, n_common, 2 * n_common]),
            ' '.join('%.3f' % T_df(p, q, x)[0] for x in [n_common // 4, n_common // 2, n_common, 2 * n_common])))
        print('     curve at n (shells %d..%d): %s' % (band[0], band[1], ' '.join('%.3f' % v for v in fc[band[0]:band[1] + 1])))
