#!/usr/bin/env python3
"""where the non-arm part of a two-map difference sits: the top 5 % of |B - A| (10 A views, standardised as in arm_readout) outside
arm_readout's closed-arm / open-arm masks, as |B - A|-weighted centroids of a 2-means split of those voxels, in map voxels (x, y, z,
0-based, box 160) and in A from the box centre, with each blob's share, extent, sign (B - A) and distance to the closed-arm centroid.
usage: diff_centroid.py A_VIEW B_VIEW"""
import sys, numpy as np
from pathlib import Path
SC = Path(__file__).resolve().parent; sys.path.insert(0, str(SC)); from model_class import read_mrc
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160'); CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
smpd, LP = 1.464, 10.
cons0 = read_mrc(str(R160 / 'toy_consensus.mrc')); n = cons0.shape[0]; f = np.fft.fftfreq(n, d=smpd)
kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G))
cons = lp(cons0); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(read_mrc(str(CL)))); d = closed - cons
thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99); mC = d > thr; mO = d < -thr
A, B = std(lp(read_mrc(sys.argv[1]))), std(lp(read_mrc(sys.argv[2]))); D = B - A
reg = supp | mC | mO; top = reg & (np.abs(D) >= np.percentile(np.abs(D[reg]), 95)); na = top & ~mC & ~mO
zz, yy, xx = np.nonzero(na); w = np.abs(D[na]); P = np.stack([xx, yy, zz], 1).astype(float)
cz, cy, cx = np.nonzero(mC); cC = np.array([cx.mean(), cy.mean(), cz.mean()])
c = P[[np.argmax(w), np.argmax(((P - P[np.argmax(w)]) ** 2).sum(1))]]
for _ in range(50):
    lab = ((P[:, None, :] - c[None]) ** 2).sum(2).argmin(1)
    c = np.array([np.average(P[lab == k], axis=0, weights=w[lab == k]) for k in range(2)])
ctr = (n - 1) / 2.
print('non-arm top-5 %% voxels: %d of %d top voxels; closed-arm centroid at voxel (%.0f, %.0f, %.0f) = (%+.0f, %+.0f, %+.0f) A from the box centre'
      % (na.sum(), top.sum(), *cC, *((cC - ctr) * smpd)))
for k in range(2):
    m = lab == k; sg = np.average(np.sign(D[na][m]), weights=w[m]); ext = np.sqrt(np.average(((P[m] - c[k]) ** 2).sum(1), weights=w[m])) * smpd
    print('blob %d: %.0f %% of the non-arm weight, centroid voxel (%.0f, %.0f, %.0f) = (%+.0f, %+.0f, %+.0f) A from the centre, rms extent %.0f A, '
          'sign of B - A %+.2f, %.0f A from the closed-arm centroid' % (k + 1, 100 * w[m].sum() / w.sum(), *c[k], *((c[k] - ctr) * smpd), ext, sg,
          np.sqrt(((c[k] - cC) ** 2).sum()) * smpd))
