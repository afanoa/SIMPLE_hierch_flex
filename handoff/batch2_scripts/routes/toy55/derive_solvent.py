#!/usr/bin/env python3
"""Solvent mask of a node from its own cells: every cell low-passed to LP, thresholded at FRAC of its maximum, the union
dilated by DIL A with a soft edge of EDGE A (so every conformation present in the node keeps its support).
usage: derive_solvent.py CELLSDIR OUT.mrc [--lp 15] [--frac 0.15] [--dil 10] [--edge 8] [--smpd 1.464]"""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import binary_dilation, distance_transform_edt
ap = argparse.ArgumentParser(); ap.add_argument('cells'); ap.add_argument('out'); ap.add_argument('--lp', type=float, default=15.); ap.add_argument('--frac', type=float, default=0.15)
ap.add_argument('--dil', type=float, default=10.); ap.add_argument('--edge', type=float, default=8.); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
fs = sorted(glob.glob(a.cells + '/cell_*.mrc')); assert fs, 'no cells in ' + a.cells
n = read_mrc(fs[0]).shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (a.lp / 1.2)**2)
union = np.zeros((n, n, n), bool)
for fn in fs:
    v = np.real(np.fft.ifftn(np.fft.fftn(read_mrc(fn).astype(np.float32)) * G)); union |= v > a.frac * v.max()
r = int(round(a.dil / a.smpd)); x = np.arange(-r, r + 1); ball = (x[:, None, None]**2 + x[None, :, None]**2 + x[None, None, :]**2) <= r * r
hard = binary_dilation(union, structure=ball); dist = distance_transform_edt(~hard) * a.smpd; mask = np.clip(1. - dist / a.edge, 0., 1.).astype(np.float32)
write_mrc(a.out, mask, a.smpd)
print('solvent mask %s: union of %d cells covers %.1f %% of the box, soft mask %.1f %% above 0.5' % (a.out, len(fs), 100. * union.mean(), 100. * (mask > 0.5).mean()))
