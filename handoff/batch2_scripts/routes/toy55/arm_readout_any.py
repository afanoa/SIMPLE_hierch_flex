#!/usr/bin/env python3
"""arm_readout.py's closed-arm / open-arm occupancies for maps of any box at the same physical size (e.g. box-64 stage-A maps): each map is
Fourier-padded or cropped to the consensus box first; a *_even_unfil.mrc argument is averaged with its *_odd_unfil.mrc.
usage: arm_readout_any.py CONSENSUS CLOSEDMAP MAP [MAP ...]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
cons_f, closed_f = sys.argv[1:3]; files = sys.argv[3:]; smpd = 1.464; LP = 10.
n = read_mrc(cons_f).shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G))
def resize(v):
    N = v.shape[0]
    if N == n: return v
    F = np.fft.fftshift(np.fft.fftn(v.astype(np.float64)))
    if N > n: lo = N // 2 - n // 2; H = F[lo:lo + n, lo:lo + n, lo:lo + n]
    else: H = np.zeros((n, n, n), complex); lo = n // 2 - N // 2; H[lo:lo + N, lo:lo + N, lo:lo + N] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(H)))
cons = lp(read_mrc(cons_f)); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(read_mrc(closed_f))); d = closed - cons; thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99)
mC = d > thr; mO = d < -thr; refC = closed[mC].mean(); refO = cons[mO].mean()
for fn in files:
    v = read_mrc(fn)
    if fn.endswith('_even_unfil.mrc'): v = 0.5 * (v + read_mrc(fn.replace('_even_unfil.mrc', '_odd_unfil.mrc')))
    v = std(lp(resize(v))); c = v[mC].mean() / refC; o = v[mO].mean() / refO
    print('%-44s closed-arm %.2f open-arm %.2f' % (os.path.basename(fn).replace('_even_unfil.mrc', ''), c, o))
