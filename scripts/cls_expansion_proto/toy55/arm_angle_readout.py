#!/usr/bin/env python3
"""Per-arm readout of a synthetic 'arms' chain: for every leaf state, the std of each arm's angles (theta, phi) over the
state's particles relative to the std over all particles; a state that localises arm k has a small ratio for arm k.
Also per region board: the arm whose angle correlates best with the cell index (the board's coordinate).
usage: arm_angle_readout.py RUNDIR GT_PARTICLES TRUTH_ARMS"""
import sys, glob, re, numpy as np
run, gtf, trf = sys.argv[1:4]
conf = {}; [conf.__setitem__(int(w[0]), int(w[1])) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')]
ang = {}; [ang.__setitem__(int(w[0]), np.array([float(x) for x in w[1:]])) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')]
narm = len(next(iter(ang.values()))) // 2; N = max(conf) + 1; A = np.full((N, 2 * narm), np.nan)
for p, c in conf.items():
    if c in ang: A[p] = ang[c]
allstd = np.nanstd(A, 0); print('arms %d, angle std over all particles (theta/phi per arm): %s' % (narm, ' '.join('%.1f' % v for v in allstd)))
for f in sorted(glob.glob(run + '/fam?/region*/level1/states_clusters.txt')):
    tag = re.search(r'(fam\d/region[^/]+)', f).group(1); st = np.loadtxt(f, dtype=int); print('== %s' % tag)
    for s in range(1, st.max() + 1):
        rows = np.where(st == s)[0] + 1; rows = rows[rows < N]; a = A[rows]; a = a[~np.isnan(a[:, 0])]
        if len(a) < 50: continue
        ratio = np.nanstd(a, 0) / allstd; k = int(np.argmin(ratio)); print('  state %d (%d particles): std ratio per angle %s  -> localises arm %d %s (ratio %.2f, mean %.1f deg)' % (s, len(a), ' '.join('%.2f' % v for v in ratio), k // 2 + 1, 'theta' if k % 2 == 0 else 'phi', ratio[k], np.mean(a[:, k])))
for b in sorted(glob.glob(run + '/fam?/region*/board')) + sorted(glob.glob(run + '/fam?/klem_fb')):
    try: L = np.loadtxt(b + '/klem_latent.txt', comments='#')
    except Exception: continue
    pind = L[:, 0].astype(int); z = L[:, 1:3]; a = A[pind]; m = ~np.isnan(a[:, 0]); tag = b.split(run.rstrip('/') + '/')[-1]
    cc = [max(abs(np.corrcoef(z[m, j], a[m, k])[0, 1]) for j in range(2)) for k in range(2 * narm)]; k = int(np.argmax(cc))
    print('board %-32s best |corr(lattice coordinate, angle)|: arm %d %s %.2f   (all: %s)' % (tag, k // 2 + 1, 'theta' if k % 2 == 0 else 'phi', cc[k], ' '.join('%.2f' % v for v in cc)))
