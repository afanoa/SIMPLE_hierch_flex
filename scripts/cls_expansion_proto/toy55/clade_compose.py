#!/usr/bin/env python3
"""Children of a data-derived tree node against the reference hierarchy: each child's particle share in each reference
clade (average-linkage clades of the true conformer maps cut at the given heights). usage: clade_compose.py CHILDREN GT TRUTH_NPZ [CUTS]"""
import sys, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
lab = np.loadtxt(sys.argv[1]).astype(int); conf = np.loadtxt(sys.argv[2], comments='#')[:, 1].astype(int)[:len(lab)]; z = np.load(sys.argv[3])
cuts = [float(c) for c in (sys.argv[4] if len(sys.argv) > 4 else '0.45,0.40').split(',')]; Z = linkage(squareform(z['D'], checks=False), 'average'); ids = list(z['ids'])
for cut in cuts:
    cl = fcluster(Z, t=cut, criterion='distance'); cmap = {k: c for k, c in zip(ids, cl)}; nc = cl.max()
    tot = np.array([sum(z['pop'][i] for i, k in enumerate(ids) if cmap[k] == c) for c in range(1, nc + 1)], float)
    print('  reference clades at %.2f: %d (population shares %s)' % (cut, nc, ' '.join('%.2f' % v for v in tot / tot.sum())))
    for c in np.unique(lab[lab > 0]):
        m = (lab == c) & (conf > 0); cc = np.bincount([cmap.get(k, 0) for k in conf[m]], minlength=nc + 1)[1:]
        print('    child %d (%d particles): clade shares %s' % (c, m.sum(), ' '.join('%.2f' % v for v in cc / max(cc.sum(), 1))))
