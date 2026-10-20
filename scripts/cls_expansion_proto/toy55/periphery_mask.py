#!/usr/bin/env python3
"""Data-derived periphery mask for the region step: the consensus map smoothed, envelope = above ENV x max, rigid core =
above CORE x max (a rigid body keeps full density in the consensus, moving parts average down), periphery = envelope
minus the dilated core, soft edge. --check DIR (arm sector masks, synthetic truth): share of each arm zone (sector beyond
RARM A) and of the inner core (< RCORE A) inside the periphery - readout only.
usage: periphery_mask.py CONS.mrc OUT.mrc [--sigma 3 --env 0.05 --core 0.40 --dilate 2 --check DIR --smpd 1.464]"""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--sigma', type=float, default=3.); ap.add_argument('--env', type=float, default=0.05)
ap.add_argument('--core', type=float, default=0.40); ap.add_argument('--dilate', type=int, default=2); ap.add_argument('--check', default=''); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--rarm', type=float, default=55.); ap.add_argument('--rcore', type=float, default=30.); a = ap.parse_args()
V = read_mrc(a.cons); s = gaussian_filter(V, a.sigma); s /= s.max(); env = s > a.env; core = binary_dilation(s > a.core, iterations=a.dilate)
P = gaussian_filter((env & ~core).astype(np.float32), 1.5); write_mrc(a.out, P.astype(np.float32), a.smpd); v3 = a.smpd ** 3 / 1000.
print('periphery mask: envelope %.0f nm^3, rigid core %.0f nm^3 (dilated), periphery %.0f nm^3 -> %s' % (env.sum() * v3, core.sum() * v3, (env & ~core).sum() * v3, a.out))
if a.check:
    B = V.shape[0]; c = B / 2; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) * a.smpd; Pb = P > 0.5
    ms = sorted(glob.glob(a.check + '/arm_*.mrc')); share = [(Pb & (read_mrc(m) > 0.5) & (r > a.rarm)).sum() / max(((read_mrc(m) > 0.5) & (r > a.rarm) & env).sum(), 1) for m in ms]
    print('   truth check: share of each arm zone (beyond %.0f A, inside the envelope) in the periphery %s; inner core (< %.0f A) in the periphery %.2f' % (a.rarm, ' '.join('%.2f' % v for v in share), a.rcore, (Pb & (r < a.rcore)).sum() / max((r < a.rcore).sum(), 1)))
