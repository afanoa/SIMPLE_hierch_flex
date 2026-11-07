#!/usr/bin/env python3
"""simsites truth maps: a RIGID RNase P core (docked conformer 49, open arm) plus three spatially DISJOINT sites, every combination
of their states one map (2 x 2 x 2 = 8 'conformers' for the simulation; make_oris_sites.py draws the sites independently per particle):
  A  deleted domain: a compact lobe of the core (ball inside the conformer-49 density holding MASSA of its mass, placed as far as
     possible from sites B and C); present (state 1) or deleted (state 0). maps160/variants/domain_mask.mrc is NOT used: it is the
     mobile arm's swept volume (86 % of the closed-arm region, 49 % of the open-arm region), i.e. the same place as site C.
  B  added blob: the blob of variants/conf62_add46 (conf62 - conf46 x mask^2, exactly the added Gaussian); present or absent.
  C  arm: open (conformer 49) or closed (docked conformer 18 inside the arm-site mask: the arm density present in one map only,
     dilated DIL px with a FALL px Gaussian falloff; conformer 49 everywhere else).
Combination id = 1 + (1 - A) + 2 B + 4 C (id 1 = domain present, no blob, arm open). Writes OUT/conf0k_sites.mrc (k = 1..8),
OUT/mask_{A,B,C}.mrc (soft site masks), OUT/consensus_sites.mrc (frequency-weighted mean), OUT/combos.txt (id A B C frequency),
OUT/snr_scale.txt (per map: projection power at the simulation mask over the toy consensus at mskdiam 200, the RNase P calibration),
CSV for merge55.py. usage: make_sites_maps.py TOY55_MAPS160 OUT CSV [--pa 0.7 --pb 0.12 --pc 0.08 --massa 0.06 --mskdiam 232]"""
import os, argparse, numpy as np
from scipy import ndimage
ap = argparse.ArgumentParser(); ap.add_argument('maps'); ap.add_argument('out'); ap.add_argument('csv')
ap.add_argument('--pa', type=float, default=0.7); ap.add_argument('--pb', type=float, default=0.12); ap.add_argument('--pc', type=float, default=0.08)
ap.add_argument('--massa', type=float, default=0.06); ap.add_argument('--mskdiam', type=float, default=232.); ap.add_argument('--dil', type=int, default=3); ap.add_argument('--fall', type=float, default=3.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--seed', type=int, default=0); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd):
    x = np.ascontiguousarray(x, np.float32); h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444; hf[54] = x.std()
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.tobytes())
s = a.smpd; M = a.maps; cons = rd(M + '/toy_consensus.mrc'); c18 = rd(M + '/docked/conf18_emd74209_m160.mrc'); c49 = rd(M + '/docked/conf49_emd74133_m160.mrc')
c46 = rd(M + '/docked/conf46_emd74259_m160.mrc'); v62 = rd(M + '/variants/conf62_add46_m160.mrc')
n = c49.shape[0]; zz, yy, xx = np.mgrid[:n, :n, :n]; r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2)
f = np.fft.fftfreq(n, d=s); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=s), indexing='ij'); k2 = kx ** 2 + ky ** 2 + kz ** 2
lpg = lambda v, res: np.fft.irfftn(np.fft.rfftn(v) * np.exp(-0.5 * k2 * (res / 1.2) ** 2), s=v.shape)
soft = lambda diam: 0.5 * (1 + np.cos(np.pi * np.clip((r - (diam / 2 / s - 4)) / 4., 0, 1)))
L49 = lpg(c49, 10); s49 = L49 > 0.2 * L49.max(); s49w = L49 > 0.1 * L49.max(); tot = c49[s49].sum()
# ---- site B: the added blob of conf62 (make_variants.py: V62 = (conf46 x msk + blob) x msk, msk = soft sphere 200 A)
blob = np.clip(v62 - c46 * soft(200.) ** 2, 0, None); bc = np.unravel_index(np.argmax(blob), blob.shape); db = np.sqrt((zz - bc[0]) ** 2 + (yy - bc[1]) ** 2 + (xx - bc[2]) ** 2) * s
print('site B blob: peak %.4f at %s, mass %.3f of conf49, residual beyond 25 A of the peak %.2e; conf49 lp10 density at the peak %.2f of its max' % (blob.max(), bc, blob.sum() / tot, blob[db > 25].max(), L49[bc] / L49.max()))
blob[db > 30] = 0.; MB = np.exp(-ndimage.distance_transform_edt(~ndimage.binary_dilation(blob > 0.05 * blob.max(), iterations=2)) ** 2 / (2 * a.fall ** 2))
# ---- site C: arm density present in one conformer only (lp 8), cleaned, dilated, Gaussian falloff
L18_8, L49_8 = lpg(c18, 8), lpg(c49, 8); mx = max(L18_8.max(), L49_8.max())
def clean(b, mn=150):
    lab, nl = ndimage.label(b); sz = np.bincount(lab.ravel()); keep = np.where(sz >= mn)[0]; keep = keep[keep > 0]; return np.isin(lab, keep)
coreC = clean((L18_8 > 0.2 * mx) & (L49_8 < 0.1 * mx)) | clean((L49_8 > 0.2 * mx) & (L18_8 < 0.1 * mx))
MC = np.exp(-ndimage.distance_transform_edt(~ndimage.binary_dilation(coreC, iterations=a.dil)) ** 2 / (2 * a.fall ** 2))
# ---- site A: a ball inside the conformer-49 density holding MASSA of its mass, as far as possible from B and C
occ = (MC > 0.05) | (MB > 0.05); far = ndimage.distance_transform_edt(~occ) * s; rng = np.random.default_rng(a.seed)
cand = np.argwhere(s49 & (far > 25.)); best = None; w49 = np.where(s49, c49, 0.)
def ball_mask(cpt, rad):
    d = np.sqrt((zz - cpt[0]) ** 2 + (yy - cpt[1]) ** 2 + (xx - cpt[2]) ** 2) * s
    return np.clip(ndimage.gaussian_filter(ndimage.binary_dilation((d <= rad) & s49w, iterations=2).astype(float), 2.) * 1.5, 0, 1)
for cpt in cand[rng.choice(len(cand), min(300, len(cand)), replace=False)]:   # screen: hard ball, mass share from the sorted radii
    d = np.sqrt((zz - cpt[0]) ** 2 + (yy - cpt[1]) ** 2 + (xx - cpt[2]) ** 2) * s; o = np.argsort(d[s49]); cm = np.cumsum(w49[s49][o]) / tot
    rad = d[s49][o][min(np.searchsorted(cm, 0.8 * a.massa), len(o) - 1)]; gap = far[(d <= rad + 4 * s) & s49w].min()
    if best is None or gap > best[0]: best = (gap, rad, tuple(int(v) for v in cpt))
cptA = best[2]; lo, hi = 4., 40.
for _ in range(12):   # bisection on the radius for the soft mask's mass share (best seed only)
    mid = 0.5 * (lo + hi); m = ball_mask(cptA, mid); lo, hi = (mid, hi) if (c49 * m)[s49].sum() / tot < a.massa else (lo, mid)
radA = hi; MA = ball_mask(cptA, radA); gapA = far[MA > 0.05].min()
print('site A deleted domain: centre %s, ball radius %.1f A, mass %.3f of conf49, soft mask > 0.5 %d voxels, min gap (mask > 0.05) to sites B/C %.1f A' % (cptA, radA, (c49 * MA)[s49].sum() / tot, (MA > 0.5).sum(), gapA))
# ---- overlaps and the arm site's readout (arm_readout.py's regions: deposited closed vs the toy consensus)
for nm, m1, m2 in (('A-B', MA, MB), ('A-C', MA, MC), ('B-C', MB, MC)): print('mask overlap %s: sum(min) / sum(smaller) %.4f' % (nm, np.minimum(m1, m2).sum() / min(m1.sum(), m2.sum())))
cs = lpg(cons, 10); supp = cs > 0.2 * cs.max(); outside = cs < 0.02 * cs.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
csd = std(cs); cld = std(lpg(c18, 10)); dd = cld - csd; thr = 0.5 * np.percentile(np.abs(dd[supp | (cld > 0.2 * cld.max())]), 99); mCr, mOr = dd > thr, dd < -thr
def readout(V): v = std(lpg(V, 10)); return v[mCr].mean() / cld[mCr].mean(), v[mOr].mean() / csd[mOr].mean()
# ---- combination maps
pa, pb, pc = a.pa, a.pb, a.pc; rows = []; maps = {}
for C_ in (0, 1):
    for B_ in (0, 1):
        for A_ in (1, 0):
            k = 1 + (1 - A_) + 2 * B_ + 4 * C_; V = c49 * (1 - C_ * MC) + C_ * c18 * MC; V = V * (1 - (1 - A_) * MA) + B_ * blob
            fr = (pa if A_ else 1 - pa) * (pb if B_ else 1 - pb) * (pc if C_ else 1 - pc); maps[k] = V; rows.append((k, A_, B_, C_, fr))
rows.sort()
pv = lambda v: np.mean([v.sum(ax).var() for ax in range(3)]); ref = pv(cons * soft(200.)); cons_s = sum(fr * maps[k] for k, _, _, _, fr in rows)
with open(a.out + '/combos.txt', 'w') as fo, open(a.out + '/snr_scale.txt', 'w') as fs:
    fo.write('# id A_domain_present B_blob_present C_arm_closed frequency (of the non-junk particles)\n')
    for k, A_, B_, C_, fr in rows:
        V = maps[k]; wr('%s/conf%02d_sites.mrc' % (a.out, k), V, s); fo.write('%d %d %d %d %.6f\n' % (k, A_, B_, C_, fr)); sc = pv(V * soft(a.mskdiam)) / ref; fs.write('conf%02d %.4f\n' % (k, sc))
        oc, oo = readout(V); print('combination %d (A %d B %d C %d) frequency %.4f: projection power x%.3f, arm readout closed %.2f open %.2f' % (k, A_, B_, C_, fr, sc, oc, oo))
    for nm, V in (('blob', 'blob'), ('broken', 'broken')): pass
wr(a.out + '/consensus_sites.mrc', cons_s, s); wr(a.out + '/mask_A.mrc', MA, s); wr(a.out + '/mask_B.mrc', MB, s); wr(a.out + '/mask_C.mrc', MC, s)
print('masks > 0.5 (%% of the conf49 support): A %.1f, B %.1f, C %.1f; consensus projection power x%.3f' % tuple([100 * ((m > 0.5) & s49).sum() / s49.sum() for m in (MA, MB, MC)] + [pv(cons_s * soft(a.mskdiam)) / ref]))
# radial coverage of the closed arm by the arm mask and the closed composite (radial_extent.py-style, closed-arm region)
edges = np.arange(60, 125, 10.); rb = np.digitize(r * s, edges) - 1; Vc = std(lpg(maps[5], 10)); line = []
for i in range(len(edges) - 1):
    m_ = mCr & (rb == i)
    if m_.sum() >= 5: line.append('%.0f-%.0f A: n %d mask %.2f closed-map %.2f' % (edges[i], edges[i + 1], m_.sum(), MC[m_].mean(), Vc[m_].mean() / max(cld[m_].mean(), 1e-6)))
print('closed-arm region by radius: ' + '; '.join(line))
with open(a.csv, 'w') as f:
    f.write('emdb_id,conformer,images,resolution,pixel,box,file\n')
    for k, A_, B_, C_, fr in rows: f.write('SITES-%d,%d,%d,0,%.3f,%d,conf%02d_sites.mrc\n' % (k, k, int(round(fr * 1e6)), s, n, k))
print('wrote %s (8 maps, masks, consensus, combos.txt, snr_scale.txt) and %s' % (a.out, a.csv))
