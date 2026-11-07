#!/usr/bin/env python3
"""usage: simple_exec prg=print_project_field projfile=P oritype=ptcl3D | dump_ptcl.py OUT.npz [NROWS]
state, corr, dfx, e1 and e2 per project row (project order), from print_project_field's key=value lines (one line per ori, those
carrying e1= or state=), saved as OUT.npz. With NROWS the row count is checked."""
import sys, numpy as np
K = ('state', 'corr', 'dfx', 'e1', 'e2'); cols = {k: [] for k in K}
for l in sys.stdin:
    if 'e1=' not in l and 'state=' not in l: continue
    d = dict(t.split('=', 1) for t in l.split() if '=' in t)
    for k in K:
        try: cols[k].append(float(d.get(k, 'nan')))
        except ValueError: cols[k].append(np.nan)
n = len(cols['e1'])
if len(sys.argv) > 2: assert n == int(sys.argv[2]), 'dumped %d rows, expected %s' % (n, sys.argv[2])
np.savez(sys.argv[1], **{k: np.array(v) for k, v in cols.items()}); print('%d rows dumped to %s' % (n, sys.argv[1]))
