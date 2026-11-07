#!/usr/bin/env python3
"""land_leafjunk.py: whole-set map-level numbers with and without a junk leaf list: pooled conformer-class shares (land_confmix outputs of
the plain leaf half-maps), particle-weighted mean closed-site occupancy and the leaf shares by closed-site band.
usage: land_leafjunk.py CONFMIX_EVEN CONFMIX_ODD ARM_RAW UNITS_LIST JUNK(comma list of leaf numbers)"""
import sys, numpy as np
ce, co, arm, ul, junk = sys.argv[1:6]; J = set(int(x) for x in junk.split(','))
n = np.array([float(l.split()[3]) for l in open(ul) if l.strip() and not l.startswith('#')]); L = len(n); keep = np.array([i + 1 not in J for i in range(L)])
def shares(f):
    S = []
    for l in open(f):
        if l.startswith('#') or 'shares' not in l: continue
        w = l.split(); i = w.index('shares'); S.append([float(w[i + 2]), float(w[i + 4]), float(w[i + 6])])
    return np.array(S)
occ = np.array([float(l.split()[1]) for l in open(arm) if l.strip() and not l.startswith('#')])
for lab, m in (('all 64 leaves', np.ones(L, bool)), ('without the %d junk leaves' % len(J), keep)):
    out = []
    for f in (ce, co): S = shares(f); out.append((n[m, None] * S[m]).sum(0) / n[m].sum())
    oc = occ[m]; nn = n[m]
    print('%-28s %7d particles | closed / intermediate / open (even) %.3f %.3f %.3f (odd) %.3f %.3f %.3f | mean closed-site %.3f | leaves >= 0.8 %.3f, 0.45-0.8 %.3f, < 0.45 %.3f' % (
        lab, nn.sum(), *out[0], *out[1], (oc * nn).sum() / nn.sum(), nn[oc >= 0.8].sum() / nn.sum(), nn[(oc >= 0.45) & (oc < 0.8)].sum() / nn.sum(), nn[oc < 0.45].sum() / nn.sum()))
print('junk leaves: mean closed-site %.3f, open-class share (even) %.3f' % ((occ[~keep] * n[~keep]).sum() / n[~keep].sum(), (n[~keep, None] * shares(ce)[~keep]).sum(0)[2] / n[~keep].sum()))
