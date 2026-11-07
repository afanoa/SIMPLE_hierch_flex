#!/usr/bin/env python3
"""Per-half pose tables with f692e's family-refined poses: every half row takes its meta row's line from the refined table
of the family that holds it (f692e/famF/refined_ptcl3D.txt, meta order); rows in no family keep the consensus line.
usage: f692e_half_poses.py F692E_DIR SIM_PREFIX LAT_PREFIX OUT_PREFIX"""
import sys, numpy as np
E, simp, latp, outp = sys.argv[1:5]
fam = np.zeros(932496, int)
for f in (1, 2): fam[np.loadtxt('%s/expand/family%d_subset.txt' % (E, f), dtype=int) - 1] = f
for h in 'AB':
    mr = np.loadtxt('%s_%s/meta_rows.txt' % (simp, h)).astype(int) - 1; base = open('%s_%s/ptcl3D_r0.txt' % (latp, h)).read().splitlines(); out = list(base)
    for f in (1, 2):
        want = {m: i for i, m in enumerate(mr) if fam[m] == f}; got = 0
        for m, l in enumerate(open('%s/fam%d/refined_ptcl3D.txt' % (E, f))):
            i = want.get(m)
            if i is not None: out[i] = l.rstrip('\n'); got += 1
        print('half %s family %d: %d rows from the refined table' % (h, f, got), flush=True)
    open('%s_%s.txt' % (outp, h), 'w').write('\n'.join(out) + '\n'); print('half %s: %d rows, %d without a family (consensus)' % (h, len(out), int((fam[mr] == 0).sum())))
