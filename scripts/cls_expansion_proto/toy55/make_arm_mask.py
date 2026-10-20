#!/usr/bin/env python3
"""Soft mask of the arm's travel: union of the closed-arm position (deposited closed conformer 18 minus consensus > thr)
and the open-arm position (the reverse), dilated and feathered. usage: make_arm_mask.py CONSENSUS CLOSEDMAP OUT.mrc [--dil 8] [--edge 6]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import binary_dilation, distance_transform_edt
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('closed'); ap.add_argument('out'); ap.add_argument('--dil', type=float, default=8.); ap.add_argument('--edge', type=float, default=6.); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
LP = 10.; n = read_mrc(a.cons).shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G))
cons = lp(read_mrc(a.cons)); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(read_mrc(a.closed))); d = closed - cons; thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99)
hard = (d > thr) | (d < -thr); r = int(round(a.dil / a.smpd)); x = np.arange(-r, r + 1); ball = (x[:, None, None]**2 + x[None, :, None]**2 + x[None, None, :]**2) <= r * r
hard = binary_dilation(hard, structure=ball); dist = distance_transform_edt(~hard) * a.smpd; mask = np.clip(1. - dist / a.edge, 0., 1.).astype(np.float32)
write_mrc(a.out, mask, a.smpd); print('arm mask %s: %.1f %% of the box above 0.5' % (a.out, 100. * (mask > 0.5).mean()))
