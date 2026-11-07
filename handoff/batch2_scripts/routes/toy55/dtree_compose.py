#!/usr/bin/env python3
"""Truth composition of the children of a data-derived tree node: per child label, particles, junk share and the leading
conformers. usage: dtree_compose.py CHILDREN GT_PARTICLES [TOP]"""
import sys, numpy as np
lab = np.loadtxt(sys.argv[1]).astype(int); conf = np.loadtxt(sys.argv[2], comments='#')[:, 1].astype(int)[:len(lab)]; top = int(sys.argv[3]) if len(sys.argv) > 3 else 3
for c in np.unique(lab[lab > 0]):
    m = lab == c; cc = np.bincount(conf[m]); nj = cc[0] if len(cc) else 0; cc[0] = 0; o = np.argsort(-cc)[:top]; nn = max(m.sum() - nj, 1)
    print('  child %d: %d particles (junk %.2f); conformers %s; %d conformers >= 5 %%' % (c, m.sum(), nj / m.sum(), ' '.join('%d:%.2f' % (k, cc[k] / nn) for k in o if cc[k] > 0), int((cc / nn >= 0.05).sum())))
