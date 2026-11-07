#!/usr/bin/env python3
"""Equal-N reconstruction labels: every group of a states file subsampled to the same particle count (min group size, or N),
the rest deselected (0); the full labels are kept as STATES_FULL. usage: neq_subsample.py STATES STATES_FULL min|N [SEED]"""
import sys, shutil, numpy as np
f, ff, spec = sys.argv[1:4]; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
shutil.copyfile(f, ff); s = np.loadtxt(f).astype(int); c = np.bincount(s)[1:]; m = int(c[c > 0].min()) if spec == 'min' else int(spec)
g = np.random.default_rng(seed); t = np.zeros_like(s)
for k in range(1, len(c) + 1):
    r = np.where(s == k)[0]
    if len(r): t[g.choice(r, min(m, len(r)), replace=False)] = k
np.savetxt(f, t, fmt='%d'); print('equal-N reconstruction: %d particles per group' % m, np.bincount(t).tolist())
