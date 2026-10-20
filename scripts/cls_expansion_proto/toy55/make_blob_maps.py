#!/usr/bin/env python3
"""Binding with NO conformational consequence: every docked apo conformer map gets a copy with a Gaussian blob of a given
protein mass added at a fixed site (the rnpA site of the holo maps, or a given centre), nothing else changed. Writes
conf2NN_blob_m160.mrc (NN = the apo conformer number) + a site mask, and a CSV with both sets (the blobbed copies carry
the same image counts, so the mixture is 50/50 and the arm swing is present in both). The blob amplitude is the mean
density inside the molecule (protein-like), its volume follows the mass at 0.81 Da/A^3.
usage: make_blob_maps.py APO_DOCKED_DIR CSV_APO OUT_DIR CSV_OUT --kda 13 [--site X,Y,Z (voxels) | --site_mask MASK]"""
import sys, os, csv, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter
ap = argparse.ArgumentParser(); ap.add_argument('apo'); ap.add_argument('csv_in'); ap.add_argument('out'); ap.add_argument('csv_out'); ap.add_argument('--kda', type=float, default=13.)
ap.add_argument('--site', default=''); ap.add_argument('--site_mask', default=''); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = list(csv.DictReader(open(a.csv_in))); maps = {int(r['conformer']): read_mrc(glob.glob('%s/conf%02d_*_m160.mrc' % (a.apo, int(r['conformer'])))[0]) for r in rows}
n = next(iter(maps.values())).shape[0]; zz, yy, xx = np.indices((n, n, n)); g = lambda v: gaussian_filter(v, 2.0)
mean_map = np.mean(list(maps.values()), 0); Lm = g(mean_map); dens = Lm > Lm.mean() + 1.5 * Lm.std()
if a.site_mask:
    mk = read_mrc(a.site_mask) > 0.5; c = np.array([np.mean(q[mk]) for q in (zz, yy, xx)])
elif a.site:
    c = np.array([float(v) for v in a.site.split(',')])
else:
    raise SystemExit('give --site or --site_mask')
vol = a.kda * 1000. / 0.81; sig = (vol / ((2 * np.pi) ** 1.5)) ** (1. / 3.) / a.smpd   # Gaussian with integral volume = vol A^3 at amplitude 1
blob = np.exp(-0.5 * ((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2) / sig ** 2); amp = float(mean_map[dens].mean()); blob *= amp
rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * a.smpd
print('blob: %.1f kDa = %.0f A^3, Gaussian sigma %.1f A, amplitude %.3f (mean molecular density), centre %s (radius %.0f A), added mass %.1f %% of the molecule' % (a.kda, vol, sig * a.smpd, amp, np.array2string(c, precision=0), float(np.hypot(np.hypot(c[2] - n / 2, c[1] - n / 2), c[0] - n / 2) * a.smpd), 100 * blob.sum() / mean_map[dens].sum()))
site = (blob > 0.25 * amp).astype(float); site = np.clip(gaussian_filter(site, 2.0) * 1.5, 0, 1); write_mrc(a.out + '/blob_site_mask.mrc', site.astype(np.float32), a.smpd)
with open(a.csv_out, 'w', newline='') as fo:
    wr = csv.DictWriter(fo, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow(r)
    for r in rows:
        k = int(r['conformer']); fn = a.out + '/conf%d_blob_m160.mrc' % (200 + k); write_mrc(fn, (maps[k] + blob).astype(np.float32), a.smpd)
        wr.writerow(dict(r, emdb_id='BLOB-%d' % k, conformer=200 + k, file=os.path.basename(fn)))
print('wrote %d blobbed maps -> %s, CSV %s (blobbed conformers = 200 + apo number)' % (len(rows), a.out, a.csv_out))
