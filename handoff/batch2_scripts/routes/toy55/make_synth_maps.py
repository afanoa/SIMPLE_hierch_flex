#!/usr/bin/env python3
"""Synthetic benchmark maps at the RNase P lattice (box 160, 1.464 A), built from REAL density (the RNase P consensus: core = its central ellipsoid, arms/subunits = a structured piece of it, rigidly placed) so the maps carry features to 3-4 A:
  arms     (IgM-like): a textured rigid core with 4 arms on hinges, every configuration has independent random arm
                       angles (in-plane +-40 deg, out-of-plane +-30 deg): many independent continuous rigid-body swings
  tethered (PDH-like): a textured core with 6 peripheral subunits present at 50 % occupancy each, every present one
                       displaced randomly within a 12 A tether radius: compositional + positional heterogeneity
Writes OUTDIR/conf%02d_synth.mrc (ids 1..NCONF), a CSV in the toy format (conformer, images, file) and TRUTH.txt with the
configuration parameters per conformer. The maps are scaled to the RMS of the RNase P toy consensus inside its mask.
usage: make_synth_maps.py MODE OUTDIR CSV TRUTH CONS.mrc [--nconf 100 --seed 0 --box 160 --smpd 1.464]"""
import sys, os, csv, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, map_coordinates
ap = argparse.ArgumentParser(); ap.add_argument('mode'); ap.add_argument('out'); ap.add_argument('csv'); ap.add_argument('truth'); ap.add_argument('cons')
ap.add_argument('--core_r', type=float, default=56.); ap.add_argument('--arm_radii', default='26,12,12'); ap.add_argument('--sub_r', type=float, default=15.); ap.add_argument('--nconf', type=int, default=100); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--box', type=int, default=160); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed); B, s = a.box, a.smpd; c = B / 2
zz, yy, xx = np.mgrid[:B, :B, :B]; X = (xx - c) * s; Y = (yy - c) * s; Z = (zz - c) * s   # Angstrom coordinates
def ellipsoid(cx, cy, cz, rx, ry, rz, R=None):
    P = np.stack([X - cx, Y - cy, Z - cz], -1)
    if R is not None: P = P @ R
    return np.exp(-0.5 * ((P[..., 0] / rx) ** 2 + (P[..., 1] / ry) ** 2 + (P[..., 2] / rz) ** 2) ** 2 * 4)   # soft super-Gaussian
def texture(mask, nblob, rb, seed):   # internal features so the core aligns: random Gaussian blobs inside the mask
    r2 = np.random.default_rng(seed); T = np.zeros_like(mask); idx = np.argwhere(mask > 0.5)
    for i in r2.choice(len(idx), nblob, replace=False):
        z0, y0, x0 = idx[i]; T += np.exp(-0.5 * (((xx - x0) * s) ** 2 + ((yy - y0) * s) ** 2 + ((zz - z0) * s) ** 2) / rb ** 2)
    return mask * (0.35 + 1.3 * T / max(T.max(), 1e-30))
def rot(ax, ang):
    ax = np.asarray(ax, float); ax /= np.linalg.norm(ax); K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]]); t = np.radians(ang)
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
cons = read_mrc(a.cons); cm = cons > cons.mean() + 1.5 * cons.std(); target_rms = np.sqrt(np.mean(cons[cm] ** 2))
from scipy.ndimage import affine_transform
def cutout(vol, centre, radii, R=None):   # real density inside a soft ellipsoid, re-centred on the box centre
    m = ellipsoid(*centre, *radii, R=R); piece = vol * m; sh = np.array([c, c, c]) - np.array([centre[2] / s + c, centre[1] / s + c, centre[0] / s + c])
    return affine_transform(piece, np.eye(3), offset=-sh, order=1)
def place(piece, Rm, t):   # rigid rotation (about the box centre) then translation t (A) of a re-centred piece
    Ri = np.linalg.inv(Rm[::-1, ::-1].T); off = np.array([c, c, c]) - Ri @ (np.array([c, c, c]) + np.array([t[2], t[1], t[0]]) / s)
    return affine_transform(piece, Ri, offset=off, order=1)
def rot_to(u0, d):   # rotation taking unit u0 to unit d
    v = np.cross(u0, d); sn = np.linalg.norm(v); cs = float(np.dot(u0, d))
    if sn < 1e-8: return np.eye(3) if cs > 0 else rot([1, 0, 0] if abs(u0[0]) < 0.9 else [0, 1, 0], 180)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]); return np.eye(3) + K + K @ K * ((1 - cs) / sn ** 2)
# the arm/subunit piece: the densest structured spot of the real map 40-60 A from the centre, cut along its radial axis
# recentre the real map on its density centroid so the core sits at the box centre
w_ = np.maximum(cons - cons.mean(), 0); cz, cy, cx = [float((w_ * q).sum() / w_.sum()) for q in (zz, yy, xx)]
cons = affine_transform(cons, np.eye(3), offset=np.array([cz, cy, cx]) - c, order=1)
def surface_radius(u, thr_frac=0.25):   # radius (A) where the smoothed density along direction u drops below thr_frac of its max
    smx = gaussian_filter(cons, 2.0); rs = np.arange(10, 110, 1.0); prof = np.array([smx[tuple(np.clip(np.rint(c + r * u[::-1] / s).astype(int), 0, B - 1))] for r in rs])
    above = np.where(prof > thr_frac * prof.max())[0]; return float(rs[above[-1]]) if len(above) else 40.
sm = gaussian_filter(cons, 3.0); rr = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); cand = np.where((rr > 40) & (rr < 60), sm, -1e30); i0 = np.unravel_index(np.argmax(cand), cand.shape)
p0 = np.array([X[i0], Y[i0], Z[i0]]); u0 = p0 / np.linalg.norm(p0); Rp = np.linalg.svd(np.outer(u0, u0))[0]
rows = []; truth = open(a.truth, 'w')
if a.mode.startswith('arms'):
    core = cons * ellipsoid(0, 0, 0, a.core_r, a.core_r, a.core_r); arm_piece = cutout(cons, p0, tuple(float(x) for x in a.arm_radii.replace('x', ',').replace(':', ',').split(',')), R=Rp); arm_piece *= target_rms / max(np.sqrt(np.mean(arm_piece[arm_piece > 0.2 * arm_piece.max()] ** 2)), 1e-30)
    hinges = [np.array(u, float) * max(surface_radius(np.array(u, float)), 30.) for u in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0])]; L = 52.   # never inside the core
    print('arms: hinge radii (A) %s' % ' '.join('%.0f' % np.linalg.norm(h) for h in hinges))
    truth.write('# conformer then per arm (4): theta_inplane phi_outofplane (deg)\n')
    for k in range(1, a.nconf + 1):
        V = core.copy(); pars = []
        for h in hinges:
            th, ph = rng.uniform(-40, 40), rng.uniform(-30, 30); u = h / np.linalg.norm(h)
            R = rot([0, 0, 1], th) @ rot(np.cross(u, [0, 0, 1]), ph); d = R @ u
            V += place(arm_piece, rot_to(u0, d), h + d * (L / 2 - 4)); pars += [th, ph]
        V *= target_rms / np.sqrt(np.mean(V[V > 0.2 * V.max()] ** 2)); write_mrc('%s/conf%02d_synth.mrc' % (a.out, k), V.astype(np.float32), s)
        rows.append(k); truth.write('%d %s\n' % (k, ' '.join('%.2f' % p for p in pars)))
elif a.mode.startswith('tethered'):
    core = cons * ellipsoid(0, 0, 0, a.core_r, a.core_r, a.core_r); sub_piece = cutout(cons, p0, (a.sub_r, a.sub_r, a.sub_r)); sub_piece *= target_rms / max(np.sqrt(np.mean(sub_piece[sub_piece > 0.2 * sub_piece.max()] ** 2)), 1e-30)
    sites = (surface_radius(np.array([1., 0, 0])) + 30.) * np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    truth.write('# conformer then per site (6): occupancy dx dy dz (A)\n')
    for k in range(1, a.nconf + 1):
        V = core.copy(); pars = []
        for p in sites:
            occ = int(rng.random() < 0.5); dv = rng.normal(0, 12. / np.sqrt(3), 3) if occ else np.zeros(3)
            if occ: V += place(sub_piece, rot([0, 0, 1], rng.uniform(0, 360)) @ rot([1, 0, 0], rng.uniform(0, 360)), p + dv)
            pars += [occ, *dv]
        V *= target_rms / np.sqrt(np.mean(V[V > 0.2 * V.max()] ** 2)); write_mrc('%s/conf%02d_synth.mrc' % (a.out, k), V.astype(np.float32), s)
        rows.append(k); truth.write('%d %s\n' % (k, ' '.join('%.2f' % p for p in pars)))
else: raise SystemExit('mode arms|tethered')
truth.close()
with open(a.csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['emdb_id', 'conformer', 'images', 'resolution', 'pixel', 'box', 'file']); w.writeheader()
    for k in rows: w.writerow(dict(emdb_id='SYN-%03d' % k, conformer=k, images=1000, resolution=0, pixel=s, box=B, file='conf%02d_synth.mrc' % k))
mean = np.mean([read_mrc('%s/conf%02d_synth.mrc' % (a.out, k)) for k in rows[:40]], 0); write_mrc('%s/consensus_synth.mrc' % a.out, mean.astype(np.float32), s)
print('%s: %d maps in %s, csv %s, truth %s, consensus of the first 40 written' % (a.mode, len(rows), a.out, a.csv, a.truth))
