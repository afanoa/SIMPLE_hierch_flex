#!/usr/bin/env python3
"""land_winconf.py: per window of a land_real.py run, the composition in deposited-conformer classes (closed >= 0.8, intermediate 0.45-0.8,
open < 0.45 by each conformer's own closed-site readout), both folds, both halves. usage: land_winconf.py RUNDIR ARM_CONFS.txt"""
import sys, json, os, numpy as np
R = json.load(open(sys.argv[1] + '/real.json')); arm = {}
for l in open(sys.argv[2]):
    if l.strip() and not l.startswith('#'): w = l.split(); arm[os.path.basename(w[0]).split('_')[0]] = float(w[1])
rc = np.array([arm.get(nm, np.nan) for nm in R['conf_names']]); cls = np.where(rc >= 0.8, 0, np.where(rc >= 0.45, 1, 2))
print('conformers: %d closed, %d intermediate, %d open' % ((cls == 0).sum(), (cls == 1).sum(), (cls == 2).sum()))
for f in (0, 1):
    for w in range(100):
        k = '%d_%d' % (f, w)
        if k not in R['windows']: break
        o = R['windows'][k]; q = [np.array(x) for x in o['qc']]; sh = [[x[cls == c].sum() for c in range(3)] for x in q]; wc = [x @ np.nan_to_num(rc) for x in q]
        print('fold %d window %2d: n %5d lam %+.3f | closed %.2f/%.2f  intermediate %.2f/%.2f  open %.2f/%.2f (halves) | weighted closed-site %.2f/%.2f' % (f, w + 1, o['n'], o['lam'],
              sh[0][0], sh[1][0], sh[0][1], sh[1][1], sh[0][2], sh[1][2], wc[0], wc[1]))
