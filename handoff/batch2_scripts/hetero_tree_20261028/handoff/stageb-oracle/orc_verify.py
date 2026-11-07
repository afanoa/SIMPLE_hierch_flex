#!/usr/bin/env python3
"""orc_verify.py DUMP [states=FILE] [poses=FILE]: a patched project's print_project_field dump carries the intended states and poses."""
import sys
import numpy as np

fd = sys.argv[1]
opts = dict(a.split('=', 1) for a in sys.argv[2:])
v = {k: [] for k in ('e1', 'e2', 'e3', 'x', 'y', 'state')}
for ln in open(fd):
    if 'e1=' not in ln:
        continue
    d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
    for k in v:
        v[k].append(float(d.get(k, 0.)))
v = {k: np.array(x) for k, x in v.items()}
msg = []
if 'states' in opts:
    s = np.loadtxt(opts['states'], dtype=int)
    assert (v['state'].astype(int) == s).all(), 'states differ'
    msg.append('states %s' % np.bincount(s).tolist())
if 'poses' in opts:
    p = np.loadtxt(opts['poses'], ndmin=2)
    for j, k in enumerate(('e1', 'e2', 'e3', 'x', 'y')):
        dv = np.abs(v[k] - p[:, j])
        if k in ('e1', 'e3'):
            dv = np.minimum(dv, 360. - dv)
        assert dv.max() < 2e-3, ('pose differs', k, dv.max())
    msg.append('poses match')
print('%s: %d rows, %s' % (fd, len(v['e1']), ', '.join(msg)))
