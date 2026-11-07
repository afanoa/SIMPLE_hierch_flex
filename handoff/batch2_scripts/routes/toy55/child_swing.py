#!/usr/bin/env python3
"""Per-group truth readout of split-test children (synthetic swing sets): n, class shares, swing share and swing mean +- sd,
swing eta2 of the grouping. usage: child_swing.py CHILDREN_A,CHILDREN_B GT [UNITS]"""
import sys, numpy as np
ch = sum(np.loadtxt(f, dtype=int) for f in sys.argv[1].split(',')); G = np.loadtxt(sys.argv[2], comments='#'); conf = G[:, 1].astype(int)[:len(ch)]
def rep(lab, tag):
    sw = (conf >= 200) & (conf <= 230) & (lab > 0); x = conf[sw] - 200.; g = lab[sw]; mu = x.mean()
    eta = sum((g == k).sum() * (x[g == k].mean() - mu) ** 2 for k in np.unique(g)) / max(((x - mu) ** 2).sum(), 1e-9)
    print('%s: swing eta2 %.3f over %d swing particles' % (tag, eta, sw.sum()))
    for k in np.unique(lab[lab > 0]):
        m = lab == k; c = conf[m]; s = (c >= 200) & (c <= 230)
        print('   %s %d: n %d  apo %.2f holo %.2f swing %.2f junk %.2f  swing %s' % (tag, k, m.sum(), ((c > 0) & (c < 101)).mean(), ((c >= 101) & (c < 200)).mean(), s.mean(), (c == 0).mean(), '%.1f +- %.1f A' % ((c[s] - 200).mean(), (c[s] - 200).std()) if s.sum() > 50 else '-'))
rep(ch, 'child')
if len(sys.argv) > 3: rep(np.loadtxt(sys.argv[3], dtype=int), 'unit')
