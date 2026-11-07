#!/usr/bin/env python3
"""A LARGE compositional variant for the RNase P toy (the assembly-intermediate regime: a whole sub-complex missing):
conf63 = the most populous conformer with everything beyond a plane removed; the plane is normal to the density's
principal axis and placed so that --frac of the mass lies beyond it (soft 3-voxel edge). CSV_OUT = CSV_IN + the row.
usage: make_big_variant.py CSV_IN DOCKED_DIR OUT_DIR CSV_OUT [--frac 0.35] [--images 186000] [--smpd 1.464]"""
import sys, os, csv, argparse, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, write_mrc, Shells
ap = argparse.ArgumentParser(); ap.add_argument('csv_in'); ap.add_argument('docked'); ap.add_argument('out'); ap.add_argument('csv_out')
ap.add_argument('--frac', type=float, default=0.35); ap.add_argument('--images', type=int, default=186000); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=200.)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = list(csv.DictReader(open(a.csv_in))); order = sorted(rows, key=lambda r: -float(r['images'])); k1 = int(order[0]['conformer'])
V = read_mrc(glob.glob('%s/conf%02d_*_m160.mrc' % (a.docked, k1))[0]); B = V.shape[0]; S = Shells(B, a.smpd); msk = S.soft_mask(a.mskdiam / 2 / a.smpd); V = V * msk
Vl = S.lowpass(V, 15); dens = Vl > Vl.mean() + 1.5 * Vl.std(); w = np.clip(Vl, 0, None) * dens
zz, yy, xx = np.indices(V.shape); P = np.stack([zz, yy, xx], -1).astype(float); X = P[dens]; wt = w[dens]; c = (X * wt[:, None]).sum(0) / wt.sum()
C = ((X - c) * wt[:, None]).T @ (X - c) / wt.sum(); ev, evec = np.linalg.eigh(C); ax = evec[:, -1]
t = (P - c) @ ax                                   # coordinate along the principal axis of the density
tv = t[dens]; o = np.argsort(tv); cum = np.cumsum(wt[o]) / wt.sum(); t0 = tv[o][min(np.searchsorted(cum, 1 - a.frac), len(o) - 1)]
keep = 1. / (1. + np.exp((t - t0) / 1.5))        # 1 on the near side of the plane, 0 beyond, ~3-voxel edge
V63 = V * keep; removed = 1 - (V63 * dens).sum() / (V * dens).sum(); vol = (dens & (keep < 0.5)).sum() * a.smpd ** 3
print('conf63 = conformer %d minus the far side of a plane normal to its principal axis (axis extent sd %.1f A, plane %.1f A from the centroid): mass removed %.1f%% (target %.0f%%), %.0f A^3 of %.0f' %
      (k1, np.sqrt(ev[-1]) * a.smpd, t0 * a.smpd, 100 * removed, 100 * a.frac, vol, dens.sum() * a.smpd ** 3))
write_mrc(a.out + '/conf63_big%d_m160.mrc' % k1, V63, a.smpd); write_mrc(a.out + '/big_keep_mask.mrc', keep * msk, a.smpd)
with open(a.csv_out, 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow(r)
    wr.writerow(dict(emdb_id='VAR-63', conformer=63, images=a.images, resolution=3.0, pixel=0.732, box=480, file='conf63_big%d_m160.mrc' % k1))
tot = float(sum(float(r['images']) for r in rows)) + a.images; print('CSV %s: conformer 63 population %.3f of the conformer particles' % (a.csv_out, a.images / tot))
