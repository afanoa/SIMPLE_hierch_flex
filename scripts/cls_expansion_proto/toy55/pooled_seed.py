#!/usr/bin/env python3
"""Seed states from the 2D sub-class-POOLED flex latent (no truth): mean latent per pose-mean sub-class (>= min members),
population-weighted k-means with K states on the pooled points, every particle inherits its sub-class's state.
Writes states.txt (one int per particle of the project, 0 = not in a kept sub-class) and prints the truth readouts.
usage: pooled_seed.py FLEXDIR PTCL2D GT --nptcl N --k 8 --out states.txt --states .. --axis .. [--min_members 20] [--tag ..]"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pool_place import read_embedding, kmeans
from latent_windows import read_gt, readout_windows, parse_axis

ap = argparse.ArgumentParser(); ap.add_argument('flexdir'); ap.add_argument('ptcl2d'); ap.add_argument('gt'); ap.add_argument('--nptcl', type=int, required=True); ap.add_argument('--k', type=int, default=8)
ap.add_argument('--out', required=True); ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True); ap.add_argument('--min_members', type=int, default=20); ap.add_argument('--tag', default='pooled')
a = ap.parse_args(); rng = np.random.default_rng(0); axis = parse_axis(a); gt = read_gt(a.gt)
_, pinds, z, _ = read_embedding(a.flexdir + '/flex_pca_embedding.bin'); zs = (z - z.mean(0)) / np.maximum(z.std(0), 1e-12)
sc = {}
for i, l in enumerate(open(a.ptcl2d)):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); c = int(float(kv.get('class', 0)))
    if c > 0: sc[i + 1] = c
sl = np.array([sc.get(int(p), 0) for p in pinds]); m = sl > 0; su, si = np.unique(sl[m], return_inverse=True); pop = np.bincount(si); ok = pop >= a.min_members
zb = np.array([zs[m][si == j].mean(0) for j in range(len(su))])
lab_s, _, _ = kmeans(zb[ok], a.k, w=pop[ok], seeds=10, rng=rng)
state_of = dict(zip(su[ok], lab_s + 1)); st = np.array([state_of.get(s, 0) for s in sl])
states = np.zeros(a.nptcl, int); states[pinds - 1] = st; np.savetxt(a.out, states, fmt='%d')
conf = np.array([gt.get(int(p), (0, 0))[0] for p in pinds]); pur = np.array([np.mean(conf[m][si == j] == np.bincount(conf[m][si == j]).argmax()) for j in range(len(su))])
print('%s: %d sub-classes (%d with >= %d members, mean majority-conformer purity %.2f), %d particles labelled, %d states' % (a.tag, len(su), ok.sum(), a.min_members, pur[ok].mean(), (st > 0).sum(), a.k))
keep = st > 0; readout_windows(st[keep], pinds[keep], gt, a.states, axis, a.k, a.tag + ' seed')
