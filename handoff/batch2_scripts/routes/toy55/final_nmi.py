#!/usr/bin/env python3
"""Truth readout of a state labelling of both halves (one label per project row, 0 = in no state): number of states,
NMI(conformer; state), NMI(reference clade; state), junk captured by the junk-majority states and real particles lost to
them. usage: final_nmi.py NAME LAB_A,LAB_B GT_A,GT_B [TRUTH_NPZ [CUT]]"""
import sys, numpy as np
name, labs, gts = sys.argv[1:4]; truth = sys.argv[4] if len(sys.argv) > 4 else ''; cut = float(sys.argv[5]) if len(sys.argv) > 5 else 0.45; clade = None
if truth:
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    z = np.load(truth); clade = dict(zip(z['ids'], fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=cut, criterion='distance')))
def nmi(tr, lab):
    m = tr > 0; t = np.unique(tr[m], return_inverse=True)[1]; g = np.unique(lab[m], return_inverse=True)[1]; J = np.zeros((t.max() + 1, g.max() + 1)); np.add.at(J, (t, g), 1); J /= J.sum()
    pt, pg = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(pt, pg)[nz])).sum() / max(-(pt[pt > 0] * np.log(pt[pt > 0])).sum(), 1e-12)
for h, (lf, gf) in enumerate(zip(labs.split(','), gts.split(','))):
    L = np.loadtxt(lf).astype(int); c = np.loadtxt(gf, comments='#')[:, 1].astype(int)[:len(L)]; m = L > 0; L, c = L[m], c[m]
    ks = np.unique(L); junkst = [k for k in ks if np.mean(c[L == k] == 0) > 0.5]; real = c > 0
    s = '%-10s half %s: %2d states (%d junk-majority); NMI(conformer; state) %.3f' % (name, 'AB'[h], len(ks), len(junkst), nmi(c, L))
    if clade is not None: s += ', NMI(clade %.2f; state) %.3f' % (cut, nmi(np.array([clade.get(v, 0) if v > 0 else 0 for v in c]), L))
    s += '; junk inside junk states %.3f, real particles inside junk states %.3f' % (np.isin(L[~real], junkst).mean() if (~real).any() else float('nan'), np.isin(L[real], junkst).mean())
    print(s, flush=True)
