#!/usr/bin/env python3
"""Multi-motion toy maps (several INDEPENDENT motions, known product populations): a fixed base map (an RNase P conformer) plus
NF rigid textured domains hinged on the base rim, 120 deg apart in the plane of the base's two long principal axes. Domain k has two
discrete states, a hinge rotation of -DELTA_k / +DELTA_k about the plane normal; its state-2 population is P_k, independent of the other
domains, so every combination's weight is the product of the per-domain marginals. Optional --cont DELTA:NBIN adds one more domain whose
hinge angle is spread evenly over NBIN bins in [-DELTA, DELTA] (a continuous motion, uniform), independent of the discrete ones.
Domains are evaluated analytically in their own rigid frame (no resampling), so both states of a domain carry identical texture.
Writes OUT/conf<id>_mm.mrc (id = digits of the discrete states, 1|2 per domain, then the continuous bin), OUT/consensus_mm.mrc (the
population-weighted mean), OUT/domain<k>_mask.mrc (both states' support, truth readout only), OUT/factors.txt (id, per-domain state,
continuous angle, weight) and the qtoy CSV (conformer,images,kind,s1,closed_arm,map; map relative to MAPROOT, s1 = continuous angle).
usage: make_mm_maps.py BASE.mrc OUT CSV MAPREL [--factors 0.5:25/0.65:25/0.8:25 (/ or ,) --cont 0 --radii 18x13x12 --seed 0]"""
import sys, os, csv, argparse, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter
ap = argparse.ArgumentParser(); ap.add_argument('base'); ap.add_argument('out'); ap.add_argument('csv'); ap.add_argument('maprel')
ap.add_argument('--factors', default='0.5:25/0.65:25/0.8:25'); ap.add_argument('--cont', default='0'); ap.add_argument('--radii', default='18x13x12')
ap.add_argument('--gap', type=float, default=-3.); ap.add_argument('--phase', type=float, default=-1.); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--images', type=float, default=1000.); ap.add_argument('--nblob', type=int, default=80); ap.add_argument('--rblob', type=float, default=2.5); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); s = a.smpd
V0 = read_mrc(a.base).astype(np.float64); B = V0.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]
w = np.maximum(V0 - 0.1 * V0.max(), 0); cen = np.array([(w * q).sum() / w.sum() for q in (xx, yy, zz)])   # xyz voxels
P = np.stack([(xx - cen[0]).ravel(), (yy - cen[1]).ravel(), (zz - cen[2]).ravel()]) * s
ev, evec = np.linalg.eigh((P * w.ravel()) @ P.T / w.sum()); ax1, ax2, nrm = evec[:, 2], evec[:, 1], evec[:, 0]   # long, mid, normal (xyz)
sm = gaussian_filter(V0, 2.0); thr = 0.15 * sm.max()
def srad(u):
    rs = np.arange(5, 110, 0.5); pts = cen[None] + rs[:, None] * u[None] / s; i = np.clip(np.rint(pts).astype(int), 0, B - 1)
    prof = sm[i[:, 2], i[:, 1], i[:, 0]]; above = np.where(prof > thr)[0]; return float(rs[above[-1]]) if len(above) else 30.
def rot(axis, ang):
    k = axis / np.linalg.norm(axis); K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]); t = np.radians(ang)
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
fac = [tuple(float(x) for x in f.split(':')) for f in a.factors.replace('/', ',').split(',') if f]; nd = len(fac)
cd, cn = (float(a.cont.split(':')[0]), int(a.cont.split(':')[1])) if ':' in a.cont else (0., 0)
ndom = nd + (1 if cn > 0 else 0); rl, r1, r2 = (float(x) for x in a.radii.split('x'))
X = (xx - cen[0]) * s; Y = (yy - cen[1]) * s; Z = (zz - cen[2]) * s
# hinge directions: ndom rim directions evenly spaced in the long-axis plane; the phase maximises the smallest rim radius (domains stay on the rim)
phs = np.radians(np.arange(0, 360. / ndom, 2.)) if a.phase < 0 else [np.radians(a.phase)]
best = max(phs, key=lambda p0: min(srad(np.cos(p0 + 2 * np.pi * k / ndom) * ax1 + np.sin(p0 + 2 * np.pi * k / ndom) * ax2) for k in range(ndom)))
U = [np.cos(best + 2 * np.pi * k / ndom) * ax1 + np.sin(best + 2 * np.pi * k / ndom) * ax2 for k in range(ndom)]; H = [u * srad(u) for u in U]
rng = np.random.default_rng(a.seed + 7); NBL, RB = a.nblob, a.rblob
blobs = [rng.uniform(-1, 1, (NBL, 3)) * np.array([rl, r1, r2]) * 0.8 for _ in range(ndom)]   # texture blobs in each domain's frame
def domain(k, ang):
    R = rot(nrm, ang); d = R @ (H[k] / np.linalg.norm(H[k])); e1 = nrm / np.linalg.norm(nrm); e2 = np.cross(d, e1)
    c = H[k] + d * (a.gap + rl); dx, dy, dz = X - c[0], Y - c[1], Z - c[2]
    q1, q2, q3 = dx * d[0] + dy * d[1] + dz * d[2], dx * e2[0] + dy * e2[1] + dz * e2[2], dx * e1[0] + dy * e1[1] + dz * e1[2]
    body = np.exp(-0.5 * ((q1 / rl) ** 2 + (q2 / r1) ** 2 + (q3 / r2) ** 2) ** 2 * 4); T = np.zeros_like(body)
    for b in blobs[k]: T += np.exp(-0.5 * ((q1 - b[0]) ** 2 + (q2 - b[1]) ** 2 + (q3 - b[2]) ** 2) / RB ** 2)
    return body * (0.35 + 1.3 * T / T.max())
def hirms(v): return np.sqrt(np.mean(v[v > 0.2 * v.max()] ** 2))
amp = hirms(np.maximum(V0, 0)) / hirms(domain(0, 0.))
dom = [{st: amp * domain(k, (-1 if st == 1 else 1) * fac[k][1]) for st in (1, 2)} for k in range(nd)]
cang = np.linspace(-cd, cd, cn) if cn > 1 else np.zeros(cn)
cdom = [amp * domain(nd, t) for t in cang]
f = np.fft.fftfreq(B, s); fr = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def band(v, lo, hi): F = np.fft.fftn(v); m = (fr >= 1. / lo) & (fr < 1. / hi); return (np.abs(F[m]) ** 2).sum()
ref = V0 + sum(dom[k][1] for k in range(nd)) + (cdom[cn // 2] if cn else 0)
print('mm: base %s, centroid %s, rim radii (A) %s, domain radii %s, amp %.3f' % (os.path.basename(a.base), np.round(cen, 1), ' '.join('%.0f' % np.linalg.norm(h) for h in H), a.radii, amp))
for k in range(nd):
    dd = dom[k][2] - dom[k][1]
    print('domain %d: P(state 2) %.2f, swing +-%.0f deg; state difference / map power: lp10 %.4f, 40-15 A %.4f, 15-10 A %.4f; tip radius %.0f A' %
          (k + 1, fac[k][0], fac[k][1], band(dd, 1e9, 10.) / band(ref, 1e9, 10.), band(dd, 40., 15.) / band(ref, 40., 15.), band(dd, 15., 10.) / band(ref, 15., 10.), np.linalg.norm(H[k]) + a.gap + 2 * rl))
if cn: dd = cdom[-1] - cdom[0]; print('continuous domain: %d bins over +-%.0f deg; end-to-end difference lp10 %.4f, 40-15 A %.4f' % (cn, cd, band(dd, 1e9, 10.) / band(ref, 1e9, 10.), band(dd, 40., 15.) / band(ref, 40., 15.)))
for k in range(ndom):
    m = (dom[k][1] + dom[k][2]) if k < nd else sum(cdom); write_mrc('%s/domain%d_mask.mrc' % (a.out, k + 1), (m > 0.1 * m.max()).astype(np.float32), s)
rows = []; cons = np.zeros_like(V0)
for st in itertools.product((1, 2), repeat=nd):
    wd = np.prod([fac[k][0] if st[k] == 2 else 1. - fac[k][0] for k in range(nd)]); Vd = V0 + sum(dom[k][st[k]] for k in range(nd))
    for j in range(max(cn, 1)):
        cid = int(''.join(str(x) for x in st) + (str(j + 1) if cn else '')); V = Vd + (cdom[j] if cn else 0); wt = wd / max(cn, 1)
        write_mrc('%s/conf%02d_mm.mrc' % (a.out, cid), V.astype(np.float32), s); cons += wt * V; rows.append((cid, st, float(cang[j]) if cn else 0., wt))
write_mrc('%s/consensus_mm.mrc' % a.out, cons.astype(np.float32), s)
with open('%s/factors.txt' % a.out, 'w') as fo:
    fo.write('# conformer %s cont_angle weight\n' % ' '.join('f%d' % (k + 1) for k in range(nd)))
    for cid, st, t, wt in rows: fo.write('%d %s %.2f %.6f\n' % (cid, ' '.join(str(x) for x in st), t, wt))
with open(a.csv, 'w', newline='') as fo:
    wr = csv.writer(fo, lineterminator='\n'); wr.writerow(['conformer', 'images', 'kind', 's1', 'closed_arm', 'map'])
    for cid, st, t, wt in rows: wr.writerow([cid, '%.3f' % (a.images * wt), 'arm', '%.2f' % t, '0.000', '%s/conf%02d_mm.mrc' % (a.maprel, cid)])
print('mm: %d maps (%d discrete combos%s), weights min %.4f max %.4f -> %s, csv %s' % (len(rows), 2 ** nd, ' x %d continuous bins' % cn if cn else '', min(r[3] for r in rows), max(r[3] for r in rows), a.out, a.csv))
