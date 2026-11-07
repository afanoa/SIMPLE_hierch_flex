#!/usr/bin/env python3
"""Prototype of audit 6.4: are two factorial sites the same variable? Posterior Spearman (what fs_rho folds at FGATE 0.5)
beside a map-space test on the sites' pole deviation maps: for each kept site its pole half maps minus the node mean (the
particle-weighted leaf maps), low-passed and least-squares scaled, inside the union of the two sites' refined masks; the
cross-half disattenuated correlation of the two deviations (cos of the angle), each deviation's own cross-half
reproducibility, and the residual test of each after the other (Fisher z with n_eff from surrogate noise correlations).
usage: pole_angle_proto.py TREE_DIR (the 1_hetero_tree dir of a finished tree with factorial/)"""
import sys, os, re, math, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mask_disc_proto import read_mrc, fourier_filter
from scipy import ndimage

def to_box(v, n):
    """Fourier crop or pad a cubic map to box n (band-limited resampling)"""
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2; h = n // 2
    if m > n:
        G = F[c - h:c - h + n, c - h:c - h + n, c - h:c - h + n]
    else:
        G = np.zeros((n, n, n), complex); cc = n // 2; hh = m // 2
        G[cc - hh:cc - hh + m, cc - hh:cc - hh + m, cc - hh:cc - hh + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))).astype(np.float32) * (n / m)**3

T = sys.argv[1]; rng = np.random.default_rng(0)
sites = {}
for l in open(os.path.join(T, 'factorial', 'table.txt')):
    w = l.split()
    if not w or l.startswith('#'): continue
    sites[w[0]] = dict(level=w[1], labels=w[2], post=w[3], rows=w[4])
kept = {}
for l in open(os.path.join(T, 'factorial', 'sites.txt')):
    w = l.split()
    if not w or l.startswith('#'): continue
    if len(w) > 10 and w[10] == 'yes': kept[w[0]] = w[9]   # refined mask path
names = [s for s in sites if s in kept and os.path.exists(os.path.join(T, 'factorial', 'maps', 'halves', s + '_pole1_even.mrc'))]
print('tree %s: kept sites with pole half maps: %s' % (T, names))
# node mean from the leaf maps (particle weighted)
leaf = [l.split() for l in open(os.path.join(T, 'hetero_tree_leaves.txt')) if l.strip()]
mean_e = None; mean_o = None; wsum = 0.
for w in leaf:
    if w[1] == 'junk': continue
    n = float(w[2]); e, smpd = read_mrc(w[5]); o, _ = read_mrc(w[6])
    mean_e = e * n if mean_e is None else mean_e + e * n; mean_o = o * n if mean_o is None else mean_o + o * n; wsum += n
mean_e /= wsum; mean_o /= wsum
nbox = mean_e.shape[0]
lp = lambda v: fourier_filter(v, smpd, 1)
me, mo = lp(mean_e), lp(mean_o); mm = 0.5 * (me + mo)
c = nbox // 2; ax = np.arange(nbox) - c; X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); sph = X**2 + Y**2 + Z**2 <= (0.45 * nbox)**2
sup = (mm > 0.1 * mm[sph].max()) & sph
dev = {}; masks = {}
for s in names:
    pe, _ = read_mrc(os.path.join(T, 'factorial', 'maps', 'halves', s + '_pole1_even.mrc')); po, _ = read_mrc(os.path.join(T, 'factorial', 'maps', 'halves', s + '_pole1_odd.mrc'))
    pe, po = lp(to_box(pe, nbox)), lp(to_box(po, nbox))
    pe *= (pe[sup] * me[sup]).sum() / max((pe[sup]**2).sum(), 1e-30); po *= (po[sup] * mo[sup]).sum() / max((po[sup]**2).sum(), 1e-30)
    dev[s] = (pe - me, po - mo)
    if os.path.exists(kept[s]):
        m, _ = read_mrc(kept[s])
        if m.shape[0] != nbox: m = ndimage.zoom(m, nbox / m.shape[0], order=1)
        masks[s] = m > 0.5
    else:
        # refined mask pruned: the site's footprint = the top 10 % of its own deviation magnitude inside the support
        mag = np.abs(0.5 * (dev[s][0] + dev[s][1])); thr = np.quantile(mag[sup], 0.9); masks[s] = (mag > thr) & sup
    print('  %s: mask %d voxels (%.1f %% of the support)%s' % (s, masks[s].sum(), 100. * (masks[s] & sup).sum() / sup.sum(), '' if os.path.exists(kept[s]) else ' [footprint: top 10 % of |deviation|, mask file pruned]'))
# noise surrogates for n_eff: half differences of the pole maps, low-passed
def neff(M):
    rs = []
    for i, s in enumerate(names):
        for t in names[i + 1:]:
            da = (dev[s][0] - dev[s][1])[M]; db = (dev[t][0] - dev[t][1])[M]
            rs.append(np.corrcoef(da, db)[0, 1])
            rs.append(np.corrcoef(da, -db[::-1])[0, 1])
    return 1. / max(np.var(rs), 1e-12)
post = {s: np.loadtxt(sites[s]['post']) for s in names}
def fz(r, n): return math.atanh(min(max(r, -0.999999), 0.999999)) * math.sqrt(max(n - 3, 1))
def corr(a, b): return float((a * b).sum() / max(math.sqrt((a * a).sum() * (b * b).sum()), 1e-30))
print('pair | posterior Spearman (fold at >= 0.5) | map angle: cross-half cos, self r_a r_b, disattenuated cos | residual z (a|b, b|a) | n_eff')
for i, a in enumerate(names):
    for b in names[i + 1:]:
        pa, pb = post[a], post[b]; ok = (pa > -1) & (pb > -1)
        ra = np.argsort(np.argsort(pa[ok])); rb = np.argsort(np.argsort(pb[ok]))
        sp = np.corrcoef(ra, rb)[0, 1]
        M = (masks[a] | masks[b]) & sph; ne = neff(M)
        Msup = sup.copy(); nsup = neff(Msup)
        ae, ao = dev[a][0][M], dev[a][1][M]; be, bo = dev[b][0][M], dev[b][1][M]
        cab = 0.5 * (corr(ae, bo) + corr(be, ao)); raa = corr(ae, ao); rbb = corr(be, bo)
        dis = cab / max(math.sqrt(max(raa * rbb, 1e-12)), 1e-12)
        def resid(xe, xo, ye, yo):
            be_ = (xe * ye).sum() / max((ye * ye).sum(), 1e-30); bo_ = (xo * yo).sum() / max((yo * yo).sum(), 1e-30)
            return corr(xe - be_ * ye, xo - bo_ * yo)
        zab = fz(resid(ae, ao, be, bo), ne - 1); zba = fz(resid(be, bo, ae, ao), ne - 1)
        print('%s %s | %.3f %s | %.3f, %.3f %.3f, %.3f (angle %.0f deg) | %.1f %.1f | %.0f' % (a, b, sp, 'FOLD' if abs(sp) >= 0.5 else 'keep', cab, raa, rbb, dis,
              math.degrees(math.acos(min(max(dis, -1), 1))), zab, zba, ne))
        ae, ao = dev[a][0][Msup], dev[a][1][Msup]; be, bo = dev[b][0][Msup], dev[b][1][Msup]
        cab = 0.5 * (corr(ae, bo) + corr(be, ao)); raa = corr(ae, ao); rbb = corr(be, bo); dis = cab / max(math.sqrt(max(raa * rbb, 1e-12)), 1e-12)
        print('      whole support: cross-half cos %.3f, self %.3f %.3f, disattenuated %.3f (angle %.0f deg), residual z %.1f %.1f, n_eff %.0f' % (cab, raa, rbb, dis,
              math.degrees(math.acos(min(max(dis, -1), 1))), fz(resid(ae, ao, be, bo), nsup - 1), fz(resid(be, bo, ae, ao), nsup - 1), nsup))
