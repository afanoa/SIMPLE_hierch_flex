#!/usr/bin/env python3
"""land_peek.py: print the window compositions (basis weights) and true window composition of a land_e1 run. usage: land_peek.py RUNDIR"""
import sys, numpy as np
Z = np.load(sys.argv[1] + '/e1_particles.npz'); Q, tb, wl, t = Z['Q'], Z['tb'], Z['wlab'], Z['t']
print('basis t: ' + ' '.join('%5.0f' % x for x in tb))
for w in range(Q.shape[0]):
    m = (wl == w) & (t >= 0); h = np.histogram(t[m], bins=np.r_[tb - 1.5, tb[-1] + 1.5])[0] / max(m.sum(), 1)
    print('w%02d unmixed ' % (w + 1) + ' '.join('%5.2f' % x for x in Q[w]) + ' | true ' + ' '.join('%5.2f' % x for x in h))
