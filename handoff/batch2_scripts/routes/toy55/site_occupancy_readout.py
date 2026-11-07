#!/usr/bin/env python3
"""Per-site readout of a synthetic 'tethered' chain: for every leaf state the mean occupancy of each site over its
particles (all particles: ~0.5); a state that selects site k has occupancy near 1 or 0 for k. Per board the site whose
occupancy correlates best with the lattice coordinate. usage: site_occupancy_readout.py RUNDIR GT_PARTICLES TRUTH_TETHERED"""
import sys, glob, re, numpy as np
run, gtf, trf = sys.argv[1:4]
conf = {}; [conf.__setitem__(int(w[0]), int(w[1])) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')]
occ = {}; [occ.__setitem__(int(w[0]), np.array([float(x) for x in w[1:]][0::4])) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')]
nsite = len(next(iter(occ.values()))); N = max(conf) + 1; O = np.full((N, nsite), np.nan)
for p, c in conf.items():
    if c in occ: O[p] = occ[c]
print('sites %d, mean occupancy over all particles: %s' % (nsite, ' '.join('%.2f' % v for v in np.nanmean(O, 0))))
for f in sorted(glob.glob(run + '/fam?/region*/level1/states_clusters.txt')):
    tag = re.search(r'(fam\d/region[^/]+)', f).group(1); st = np.loadtxt(f, dtype=int); print('== %s' % tag)
    for s in range(1, st.max() + 1):
        rows = np.where(st == s)[0] + 1; rows = rows[rows < N]; o = O[rows]; o = o[~np.isnan(o[:, 0])]
        if len(o) < 50: continue
        m = o.mean(0); k = int(np.argmax(np.abs(m - 0.5))); print('  state %d (%d particles): site occupancies %s -> site %d %s (%.2f)' % (s, len(o), ' '.join('%.2f' % v for v in m), k + 1, 'present' if m[k] > 0.5 else 'absent', m[k]))
for b in sorted(glob.glob(run + '/fam?/region*/board')) + sorted(glob.glob(run + '/fam?/klem_fb')) + [run + '/parent']:
    try: L = np.loadtxt(b + '/klem_latent.txt', comments='#')
    except Exception: continue
    pind = L[:, 0].astype(int); z = L[:, 1:3]; o = O[pind]; m = ~np.isnan(o[:, 0]); tag = b.split(run.rstrip('/') + '/')[-1]
    cc = [max(abs(np.corrcoef(z[m, j], o[m, k])[0, 1]) for j in range(2)) for k in range(nsite)]; k = int(np.argmax(cc))
    print('board %-32s best |corr(lattice coordinate, site occupancy)|: site %d %.2f   (all: %s)' % (tag, k + 1, cc[k], ' '.join('%.2f' % v for v in cc)))
