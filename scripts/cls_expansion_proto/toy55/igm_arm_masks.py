#!/usr/bin/env python3
"""Per-arm region masks of an IgM-like set from make_igm_maps' hinge layout (same seed -> same hinge azimuths): arm k =
the azimuthal sector of width 360/NARM centred on its hinge, beyond RMIN A from the axis, soft edges; plus all_arms.mrc.
--pairs: two Fabs per hinge -> arm_{2k-1} and arm_{2k} are both written as hinge k's sector (the truth numbering).
usage: igm_arm_masks.py OUTDIR [--narm 5 --seed 0 --rmin 35 --box 160 --smpd 1.464 --pairs]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import write_mrc
from scipy.ndimage import gaussian_filter
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--narm', type=int, default=5); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--rmin', type=float, default=35.); ap.add_argument('--box', type=int, default=160); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--pairs', action='store_true'); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); B, s = a.box, a.smpd; c = B / 2; zz, yy, xx = np.mgrid[:B, :B, :B]; X, Y, Z = (xx - c) * s, (yy - c) * s, (zz - c) * s
ph0 = np.radians(np.random.default_rng(a.seed).uniform(0, 360 / a.narm))   # the first draw of make_igm_maps' generator
az = np.arctan2(Y, X); r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); half = np.pi / a.narm; allm = np.zeros((B, B, B), np.float32)
for k in range(a.narm):
    d = np.angle(np.exp(1j * (az - ph0 - 2 * np.pi * k / a.narm))); m = ((np.abs(d) < half) & (r > a.rmin) & (r < c * s - 2)).astype(np.float32)
    m = gaussian_filter(m, 2.0); allm = np.maximum(allm, m)
    for j in ((2 * k + 1, 2 * k + 2) if a.pairs else (k + 1,)): write_mrc('%s/arm_%d.mrc' % (a.out, j), m, s)
write_mrc('%s/all_arms.mrc' % a.out, allm, s); print('%d arm sector masks + all_arms.mrc in %s (hinge azimuth 0 = %.1f deg)' % (a.narm, a.out, np.degrees(ph0)))
