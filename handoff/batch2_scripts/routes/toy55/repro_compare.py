#!/usr/bin/env python3
"""Board reproduction check: compare klem_latent.txt (pind z1 z2 argmax maxr gridmass) and the cell maps of a test
board against a reference board. usage: repro_compare.py REF_DIR TEST_DIR"""
import sys, glob, os, numpy as np
ref, tst = sys.argv[1:3]
A = np.loadtxt(ref + '/klem_latent.txt', comments='#'); B = np.loadtxt(tst + '/klem_latent.txt', comments='#')
assert A.shape == B.shape and (A[:, 0] == B[:, 0]).all(), 'particle lists differ'
print('latent: %d particles; identical rows %.4f; argmax identical %.4f; max |dz| %.2e; max |d maxr| %.2e; max |d gridmass| %.2e'
      % (len(A), np.mean((A == B).all(1)), np.mean(A[:, 3] == B[:, 3]), np.abs(A[:, 1:3] - B[:, 1:3]).max(), np.abs(A[:, 4] - B[:, 4]).max(), np.abs(A[:, 5] - B[:, 5]).max()))
def rd(f):
    h = np.fromfile(f, np.int32, 256); nx, ny, nz = h[:3]; nsym = h[23]
    return np.fromfile(f, np.float32, offset=1024 + nsym).reshape(nz, ny, nx)
cr = sorted(glob.glob(ref + '/cells/*.mrc')); worst = (0., '')
for f in cr:
    g = tst + '/cells/' + os.path.basename(f)
    if not os.path.isfile(g): print('missing cell map', os.path.basename(f)); continue
    a, b = rd(f), rd(g); d = np.abs(a - b).max() / max(a.std(), 1e-30)
    if d >= worst[0]: worst = (d, os.path.basename(f))
print('cell maps: %d compared; worst max|diff|/std %.2e (%s)' % (len(cr), worst[0], worst[1]))
print('VERDICT: %s' % ('IDENTICAL' if np.mean((A == B).all(1)) == 1.0 and worst[0] == 0 else 'DIFFERENT (see above)'))
