#!/usr/bin/env python3
"""IgM-like synthetic maps at the RNase P lattice (box 160, 1.464 A): a large asymmetric core that carries most of the
mass (a filled, internally textured ellipsoid plus the RNase P consensus density for secondary-structure detail) and
NARM Fab-sized arms (filled textured ellipsoids) hinged on the core rim in the xy plane at 360/NARM spacing; every
configuration draws independent in-plane / out-of-plane hinge angles per arm (a product space, as IgM's Fabs).
Writes OUTDIR/conf%02d_synth.mrc, consensus_synth.mrc, the toy CSV, TRUTH (per conformer: theta phi per arm) and
OUTDIR/snr_scale.txt = projection power of the consensus over that of the RNase P toy consensus: simulate_particles sets
the noise from each image's own variance, so snr x this factor keeps the RNase P absolute noise level (bigger particle,
more signal at the same dose).
--pairs: two Fabs on every hinge (2 NARM arms, IgM's Fab pairs), each Fab with its own angles; in-plane = +-SPLAY (the pair's
opening) + U(-THETA, THETA), out-of-plane U(-PHI, PHI); truth columns per Fab (hinge k -> arms 2k-1, 2k).
usage: make_igm_maps.py OUTDIR CSV TRUTH CONS.mrc [--narm 5 --core_radii 54x46x40 --arm_radii 36x22x20 --arm_offset 0.55 --theta 25 --phi 30 --nconf 200 --pairs --splay 18]"""
import sys, os, csv, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, affine_transform
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('csv'); ap.add_argument('truth'); ap.add_argument('cons')
ap.add_argument('--narm', type=int, default=5); ap.add_argument('--core_radii', default='54x46x40'); ap.add_argument('--arm_radii', default='36x22x20'); ap.add_argument('--arm_offset', type=float, default=0.55)
ap.add_argument('--theta', type=float, default=25.); ap.add_argument('--phi', type=float, default=30.); ap.add_argument('--rna_w', type=float, default=0.7); ap.add_argument('--pairs', action='store_true'); ap.add_argument('--splay', type=float, default=18.)
ap.add_argument('--nconf', type=int, default=200); ap.add_argument('--images', type=int, default=1000); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--box', type=int, default=160); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed); B, s = a.box, a.smpd; c = B / 2
zz, yy, xx = np.mgrid[:B, :B, :B]; X = (xx - c) * s; Y = (yy - c) * s; Z = (zz - c) * s
radii = lambda t: tuple(float(x) for x in t.replace(',', 'x').split('x'))
def ellipsoid(cx, cy, cz, rx, ry, rz):
    return np.exp(-0.5 * (((X - cx) / rx) ** 2 + ((Y - cy) / ry) ** 2 + ((Z - cz) / rz) ** 2) ** 2 * 4)
def texture(mask, nblob, rb, seed):   # filled body with internal Gaussian blobs (features to ~2 rb)
    r2 = np.random.default_rng(seed); T = np.zeros_like(mask); idx = np.argwhere(mask > 0.5)
    for i in r2.choice(len(idx), nblob, replace=False):
        z0, y0, x0 = idx[i]; T += np.exp(-0.5 * (((xx - x0) * s) ** 2 + ((yy - y0) * s) ** 2 + ((zz - z0) * s) ** 2) / rb ** 2)
    return mask * (0.35 + 1.3 * T / max(T.max(), 1e-30))
def hirms(v): return np.sqrt(np.mean(v[v > 0.2 * v.max()] ** 2))
def rot(ax, ang):
    ax = np.asarray(ax, float); ax /= np.linalg.norm(ax); K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]]); t = np.radians(ang)
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
def place(piece, Rm, t):   # rigid rotation about the box centre, then translation t (A, xyz)
    Ri = np.linalg.inv(Rm[::-1, ::-1].T); off = np.array([c, c, c]) - Ri @ (np.array([c, c, c]) + np.array([t[2], t[1], t[0]]) / s)
    return affine_transform(piece, Ri, offset=off, order=1)
def rot_to(u0, d):
    v = np.cross(u0, d); sn = np.linalg.norm(v); cs = float(np.dot(u0, d))
    if sn < 1e-8: return np.eye(3) if cs > 0 else rot([0, 0, 1], 180)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]); return np.eye(3) + K + K @ K * ((1 - cs) / sn ** 2)
cons = read_mrc(a.cons); w_ = np.maximum(cons - cons.mean(), 0); cz, cy, cx = [float((w_ * q).sum() / w_.sum()) for q in (zz, yy, xx)]
cons = affine_transform(cons, np.eye(3), offset=np.array([cz, cy, cx]) - c, order=1); target_rms = hirms(cons)
cr = radii(a.core_radii); body = texture(ellipsoid(0, 0, 0, *cr), 160, 5., a.seed + 11); rna = np.maximum(cons, 0) * ellipsoid(0, 0, 0, *[0.95 * r for r in cr])
core = body / hirms(body) + a.rna_w * rna / hirms(rna)
ar = radii(a.arm_radii); arm = texture(ellipsoid(0, 0, 0, *ar), 30, 4., a.seed + 12); arm /= hirms(arm)   # long axis along +x, centred
u0 = np.array([1., 0, 0]); smc = gaussian_filter(core, 2.0)
def surface_radius(u, thr=0.25):
    rs = np.arange(10, 110, 1.0); prof = np.array([smc[tuple(np.clip(np.rint(c + r * u[::-1] / s).astype(int), 0, B - 1))] for r in rs])
    above = np.where(prof > thr * prof.max())[0]; return float(rs[above[-1]]) if len(above) else 40.
ph0 = np.radians(rng.uniform(0, 360 / a.narm)); hinges = []
for k in range(a.narm):
    u = np.array([np.cos(ph0 + 2 * np.pi * k / a.narm), np.sin(ph0 + 2 * np.pi * k / a.narm), 0.]); hinges.append(u * surface_radius(u))
print('igm: core radii %s, %d arms %s, hinge radii (A) %s, arm tip <= %.0f A' % (a.core_radii, a.narm, a.arm_radii, ' '.join('%.0f' % np.linalg.norm(h) for h in hinges), max(np.linalg.norm(h) for h in hinges) + (a.arm_offset + 1) * ar[0]))
nfab = a.narm * (2 if a.pairs else 1); pc, pa = (core ** 2).sum(), (arm ** 2).sum(); print('power fractions (sum V^2, one configuration): core %.2f, each of %d arms %.3f' % (pc / (pc + nfab * pa), nfab, pa / (pc + nfab * pa)))
rows = []; truth = open(a.truth, 'w'); truth.write('# conformer then per arm (%d): theta_inplane phi_outofplane (deg)\n' % nfab)
for k in range(1, a.nconf + 1):
    V = core.copy(); pars = []
    for h, off in [(h, o) for h in hinges for o in ((-a.splay, a.splay) if a.pairs else (0.,))]:
        th, ph = off + rng.uniform(-a.theta, a.theta), rng.uniform(-a.phi, a.phi); u = h / np.linalg.norm(h)
        d = rot([0, 0, 1], th) @ rot(np.cross(u, [0, 0, 1]), ph) @ u
        V += place(arm, rot_to(u0, d), h + d * a.arm_offset * ar[0]); pars += [th, ph]
    V *= target_rms / hirms(V); write_mrc('%s/conf%02d_synth.mrc' % (a.out, k), V.astype(np.float32), s)
    rows.append(k); truth.write('%d %s\n' % (k, ' '.join('%.2f' % p for p in pars)))
truth.close()
with open(a.csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['emdb_id', 'conformer', 'images', 'resolution', 'pixel', 'box', 'file']); w.writeheader()
    for k in rows: w.writerow(dict(emdb_id='SYN-%03d' % k, conformer=k, images=a.images, resolution=0, pixel=s, box=B, file='conf%02d_synth.mrc' % k))
mean = np.mean([read_mrc('%s/conf%02d_synth.mrc' % (a.out, k)) for k in rows[:40]], 0); write_mrc('%s/consensus_synth.mrc' % a.out, mean.astype(np.float32), s)
pv = lambda v: np.mean([v.sum(ax).var() for ax in range(3)]); ref = read_mrc(a.cons); first = read_mrc('%s/conf%02d_synth.mrc' % (a.out, rows[0]))
sc = pv(first) / pv(ref); open('%s/snr_scale.txt' % a.out, 'w').write('%.4f\n' % sc)
print('igm: %d maps in %s, csv %s, truth %s; projection power vs the RNase P consensus x%.2f (snr_scale.txt)' % (len(rows), a.out, a.csv, a.truth, sc))
