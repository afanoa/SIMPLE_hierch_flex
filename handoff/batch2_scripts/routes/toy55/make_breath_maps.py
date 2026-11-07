#!/usr/bin/env python3
"""Synthetic breathing (distributed contraction) of one conformer: the periphery contracts smoothly toward the density centroid
while the core stays fixed, x -> c + (1 - alpha g(r)) (x - c), g = smoothstep from 0 at R0 to 1 at R1 (A), alpha from 0 to
AMAX in NSTEPS equal steps (step 0 = the input map). Maps by backward mapping (map_coordinates, cubic). Writes
OUT/conf<ID0+k>_breath_m160.mrc and OUT/truth_breath.txt (step, conformer, alpha, inward displacement at R1 and at the
density's 95th-percentile radius; A). usage: make_breath_maps.py MAP OUT [--amax 0.10 --nsteps 11 --r0 25 --r1 55 --id0 300]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, map_coordinates
ap = argparse.ArgumentParser(); ap.add_argument('map'); ap.add_argument('out'); ap.add_argument('--amax', type=float, default=0.10); ap.add_argument('--nsteps', type=int, default=11)
ap.add_argument('--r0', type=float, default=25.); ap.add_argument('--r1', type=float, default=55.); ap.add_argument('--id0', type=int, default=300); ap.add_argument('--smpd', type=float, default=1.464)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
V = read_mrc(a.map).astype(np.float32); B = V.shape[0]; Vl = gaussian_filter(V, 1.5); dens = Vl > Vl.mean() + 1.5 * Vl.std()
idx = np.indices(V.shape).astype(np.float32); w = np.maximum(V, 0) * dens; c = np.array([(idx[i] * w).sum() / w.sum() for i in range(3)])
X = idx - c[:, None, None, None]; r = np.sqrt((X ** 2).sum(0)) * a.smpd; s = np.clip((r - a.r0) / (a.r1 - a.r0), 0, 1); g = s * s * (3 - 2 * s)
r95 = np.percentile(r[dens], 95); g95 = float(np.interp(r95, [a.r0, a.r1], [0, 1]) if r95 < a.r1 else 1.); lines = []
for k in range(a.nsteps):
    al = a.amax * k / (a.nsteps - 1)
    src = c[:, None, None, None] + X / np.maximum(1 - al * g, 1e-3)   # the output point x came from the input point at radius r / (1 - alpha g)
    out = map_coordinates(V, src.reshape(3, -1), order=3, mode='constant').reshape(V.shape) if k else V
    write_mrc('%s/conf%d_breath_m160.mrc' % (a.out, a.id0 + k), out.astype(np.float32), a.smpd)
    lines.append('%d %d %.4f %.2f %.2f' % (k, a.id0 + k, al, al * a.r1, al * g95 * r95))
with open(a.out + '/truth_breath.txt', 'w') as f:
    f.write('# step conformer alpha inward_at_R1_A inward_at_r95_A (centre %s, R0 %.0f A, R1 %.0f A, r95 %.1f A)\n' % (np.round(c, 1), a.r0, a.r1, r95))
    f.write('\n'.join(lines) + '\n')
print('breathing: centre %s, r95 %.1f A; alpha 0 .. %.2f in %d steps -> inward %.1f A at R1 (%.0f A), %.1f A at r95' % (np.round(c, 1), r95, a.amax, a.nsteps, a.amax * a.r1, a.r1, a.amax * g95 * r95))
