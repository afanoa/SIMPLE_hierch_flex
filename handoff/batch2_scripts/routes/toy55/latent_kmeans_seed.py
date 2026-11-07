#!/usr/bin/env python3
"""Truth-free seed states from the per-particle flex latent: standardise, drop the 5 % largest-norm particles from the
centroid fit (junk/outliers), population-balanced k-means with K states, every embedded particle assigned to its
nearest centroid. Writes states.txt (one int per project particle, 0 = not embedded) and prints the truth readouts.
usage: latent_kmeans_seed.py FLEXDIR GT --nptcl N --k 8 --out states.txt --states .. --axis .. [--ndim 0] [--tag ..]"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pool_place import read_embedding, kmeans
from latent_windows import read_gt, readout_windows, parse_axis
ap = argparse.ArgumentParser(); ap.add_argument('flexdir'); ap.add_argument('gt'); ap.add_argument('--nptcl', type=int, required=True); ap.add_argument('--k', type=int, default=8); ap.add_argument('--out', required=True)
ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True); ap.add_argument('--ndim', type=int, default=0); ap.add_argument('--tag', default='kseed')
a = ap.parse_args(); rng = np.random.default_rng(0); axis = parse_axis(a); gt = read_gt(a.gt)
_, pinds, z, _ = read_embedding(a.flexdir + '/flex_pca_embedding.bin'); zs = (z - z.mean(0)) / np.maximum(z.std(0), 1e-12)
if a.ndim: zs = zs[:, :a.ndim]
nrm = np.linalg.norm(zs, axis=1); core = nrm <= np.percentile(nrm, 95)
sub = rng.choice(np.where(core)[0], min(60000, core.sum()), replace=False)
lab_sub, cent, _ = kmeans(zs[sub], a.k, seeds=10, rng=rng)
d2 = ((zs[:, None, :] - cent[None, :, :]) ** 2).sum(2); st = d2.argmin(1) + 1
states = np.zeros(a.nptcl, int); states[pinds - 1] = st; np.savetxt(a.out, states, fmt='%d')
print('%s: %d particles embedded, %d dims, %d states (centroids from %d core particles); sizes %s' % (a.tag, len(pinds), zs.shape[1], a.k, len(sub), np.bincount(st)[1:]))
readout_windows(st, pinds, gt, a.states, axis, a.k, a.tag + ' seed')
