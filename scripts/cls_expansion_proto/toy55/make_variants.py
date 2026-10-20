#!/usr/bin/env python3
"""Compositional variants for the RNase P toy: from the docked conformer maps find the MOBILE DOMAIN (largest connected
component of the across-conformer variance inside the density) and write two extra truth states:
  conf61 = most populous conformer with the mobile domain deleted (domain absent)
  conf62 = second most populous conformer with an extra bound density (Gaussian blob on the surface, away from the domain)
plus the CSV with the two rows added (image counts set the populations). The maps go to a separate directory so the
10-conformer simulation does not pick them up.
usage: make_variants.py CSV_IN DOCKED_DIR CONSENSUS OUT_DIR CSV_OUT [--smpd 1.464] [--images 160000,100000]"""
import sys, os, csv, argparse, glob
import numpy as np
from scipy.ndimage import gaussian_filter, label, binary_dilation, binary_erosion
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, write_mrc, Shells

ap = argparse.ArgumentParser(); ap.add_argument('csv_in'); ap.add_argument('docked'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('csv_out')
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--images', default='160000,100000'); ap.add_argument('--mskdiam', type=float, default=200.)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = list(csv.DictReader(open(a.csv_in))); order = sorted(rows, key=lambda r: -float(r['images']))
Vc = read_mrc(a.cons); B = Vc.shape[0]; S = Shells(B, a.smpd); msk = S.soft_mask(a.mskdiam / 2 / a.smpd)
maps = {int(r['conformer']): read_mrc(glob.glob('%s/conf%02d_*_m160.mrc' % (a.docked, int(r['conformer'])))[0]) * msk for r in rows}
w = np.array([float(r['images']) for r in rows]); w /= w.sum()
L = {k: S.lowpass(v, 15) for k, v in maps.items()}; mean = sum(wi * L[int(r['conformer'])] for wi, r in zip(w, rows))
var = sum(wi * (L[int(r['conformer'])] - mean) ** 2 for wi, r in zip(w, rows))
Vl = S.lowpass(Vc, 15); dens = Vl > Vl.mean() + 1.5 * Vl.std(); dens_d = binary_dilation(dens, iterations=2)
thr = np.percentile(var[dens_d], 85); hi = (var > thr) & dens_d
lab, nlab = label(hi); sizes = np.bincount(lab.ravel())[1:]; big = np.argmax(sizes) + 1; comp = lab == big
comp = binary_dilation(comp, iterations=3); dom = np.clip(gaussian_filter(comp.astype(float), 2.) * 1.5, 0, 1)
zz, yy, xx = np.indices(Vc.shape); cen = np.array([np.mean(q[comp]) for q in (zz, yy, xx)]); mol = np.array([np.mean(q[dens]) for q in (zz, yy, xx)])
vol_dom = comp.sum() * a.smpd ** 3; vol_mol = dens.sum() * a.smpd ** 3
print('mobile domain: %d voxels = %.0f A^3 = %.1f%% of the molecular volume (%.0f A^3); centroid %s, %.1f A from the molecule centroid; variance share inside %.2f' %
      (comp.sum(), vol_dom, 100 * vol_dom / vol_mol, vol_mol, np.array2string(cen, precision=1), np.linalg.norm(cen - mol) * a.smpd, var[comp].sum() / var[dens_d].sum()))
# conf61: domain deleted from the most populous conformer
k1 = int(order[0]['conformer']); V61 = maps[k1] * (1 - dom)
# conf62: extra density on the second most populous conformer, on the surface far from the domain, inside the mask
k2 = int(order[1]['conformer']); surf = dens & ~binary_erosion(dens, iterations=2)
r_box = np.hypot(np.hypot(xx - B / 2, yy - B / 2), zz - B / 2) * a.smpd
cand = surf & (r_box < a.mskdiam / 2 - 25.)
dist_dom = np.sqrt(sum((q - c) ** 2 for q, c in zip((zz, yy, xx), cen))); p = np.unravel_index(np.argmax(np.where(cand, dist_dom, -1)), Vc.shape)
p = np.array(p, float); outward = (p - mol) / max(np.linalg.norm(p - mol), 1e-6); c_blob = p + outward * (8. / a.smpd)
sig = 8. / a.smpd; blob = np.exp(-0.5 * sum((q - c) ** 2 for q, c in zip((zz, yy, xx), c_blob)) / sig ** 2)
amp = float(np.mean(maps[k2][dens])); blob *= amp; V62 = (maps[k2] + blob) * msk
print('conf61 = conformer %d minus the domain (mass removed %.1f%%); conf62 = conformer %d plus a blob at %s (%.1f A from the domain, mass added %.1f%%)' %
      (k1, 100 * (maps[k1] * dom)[dens].sum() / maps[k1][dens].sum(), k2, np.array2string(c_blob, precision=0), np.linalg.norm(c_blob - cen) * a.smpd, 100 * blob.sum() / maps[k2][dens].sum()))
write_mrc(a.out + '/conf61_del%d_m160.mrc' % k1, V61, a.smpd); write_mrc(a.out + '/conf62_add%d_m160.mrc' % k2, V62, a.smpd); write_mrc(a.out + '/domain_mask.mrc', dom, a.smpd)
n61, n62 = [int(v) for v in a.images.split(',')]
with open(a.csv_out, 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader()
    for r in rows: wr.writerow(r)
    wr.writerow(dict(emdb_id='VAR-61', conformer=61, images=n61, resolution=3.0, pixel=0.732, box=480, file='conf61_del%d_m160.mrc' % k1))
    wr.writerow(dict(emdb_id='VAR-62', conformer=62, images=n62, resolution=3.0, pixel=0.732, box=480, file='conf62_add%d_m160.mrc' % k2))
tot = float(sum(float(r['images']) for r in rows)) + n61 + n62
print('CSV %s: conformer populations 61 %.3f, 62 %.3f of the conformer particles' % (a.csv_out, n61 / tot, n62 / tot))
