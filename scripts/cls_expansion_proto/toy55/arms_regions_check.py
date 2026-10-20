#!/usr/bin/env python3
"""Where are the arms, and where are the derived regions? Arm k's travel = |map(conformer with max theta_k) - map(min theta_k)|
from the synthetic conformer maps; region masks = the chain's regions/regions_motion of each family. Prints the centroid
and the overlap (fraction of the arm travel inside each region mask). usage: arms_regions_check.py RUNDIR SYNTHMAPS TRUTH_ARMS"""
import sys, glob, os, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
run, maps, trf = sys.argv[1:4]; ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
mp = lambda c: read_mrc('%s/conf%02d_synth.mrc' % (maps, c)).astype(np.float32)
narm = len(next(iter(ang.values()))) // 2; confs = sorted(c for c in ang if os.path.isfile('%s/conf%02d_synth.mrc' % (maps, c)))
# arm k's travel = voxels whose density depends on arm k's angles: per-voxel least squares of the 100 maps on the 2*narm angles
X = np.array([np.r_[1.0, ang[c]] for c in confs]); Y = np.array([mp(c).ravel() for c in confs]); beta = np.linalg.lstsq(X, Y, rcond=None)[0]; n = int(round(Y.shape[1] ** (1 / 3)))
arms = []
from model_class import write_mrc
for k in range(narm):
    w = np.sqrt(beta[1 + 2 * k] ** 2 * X[:, 1 + 2 * k].std() ** 2 + beta[2 + 2 * k] ** 2 * X[:, 2 + 2 * k].std() ** 2).reshape(n, n, n); m = w > 0.3 * w.max()
    idx = np.argwhere(m); cen = (idx * w[m][:, None]).sum(0) / w[m].sum(); arms.append(m); write_mrc('%s/arm%d_travel.mrc' % (run, k + 1), m.astype(np.float32), 1.464)
    print('arm %d travel (voxels whose density depends on its angles): centroid %s, %d voxels -> %s/arm%d_travel.mrc' % (k + 1, np.round(cen, 1), m.sum(), run, k + 1))
for mk in sorted(glob.glob(run + '/fam?/regions*/region_*.mrc')) + sorted(glob.glob(run + '/regions_parent/region_*.mrc')):
    m = read_mrc(mk).astype(np.float32) > 0.5; idx = np.argwhere(m); cen = idx.mean(0) if len(idx) else np.full(3, np.nan)
    print('%-45s centroid %s, %6d voxels, arm travel inside: %s' % (mk.split(run.rstrip('/') + '/')[-1], np.round(cen, 1), m.sum(), ' '.join('%.2f' % ((a & m).sum() / max(a.sum(), 1)) for a in arms)))
