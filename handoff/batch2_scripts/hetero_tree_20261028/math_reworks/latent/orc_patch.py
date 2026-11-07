#!/usr/bin/env python3
"""orc_patch.py IN OUT [states=FILE] [poses=FILE] [check=DUMP]: a copy of a SIMPLE project with the ptcl3D state and/or pose
slots replaced and every other byte kept (simple_binoris: a header of 20 segments x 5 int64 [fromto(2), bytes per record,
records, first data byte (1-based)]; ptcl3D is segment 6, one record of 64 float32 per particle, slots of simple_defs_ori).
check=DUMP compares the unpatched record with that print_project_field dump first."""
import sys
import numpy as np

I_E1, I_E2, I_E3, I_STATE, I_X, I_Y = 8, 9, 10, 23, 27, 30
fin, fout = sys.argv[1:3]
opts = dict(a.split('=', 1) for a in sys.argv[3:])
buf = bytearray(open(fin, 'rb').read())
hdr = np.frombuffer(bytes(buf[:800]), '<i8').reshape(20, 5)
f1, f2, nb, nr, first = (int(v) for v in hdr[5])
assert nb == 256 and f1 == 1 and f2 == nr and first > 800, hdr[5]
off = first - 1
rec = np.frombuffer(bytes(buf[off:off + nr * nb]), '<f4').reshape(nr, 64).copy()
if 'check' in opts:
    vals = {k: [] for k in ('e1', 'e2', 'e3', 'x', 'y', 'state')}
    for ln in open(opts['check']):
        if 'e1=' not in ln:
            continue
        dd = dict(t.split('=', 1) for t in ln.split() if '=' in t)
        for k in vals:
            vals[k].append(float(dd.get(k, 0.)))
    for k, s in zip(('e1', 'e2', 'e3', 'x', 'y', 'state'), (I_E1, I_E2, I_E3, I_X, I_Y, I_STATE)):
        dv = np.abs(rec[:, s - 1] - np.array(vals[k])).max()
        assert dv < 1e-3, ('slot check failed', k, dv)
    print('slot layout matches the dump over %d records' % nr)
if 'states' in opts:
    s = np.loadtxt(opts['states'], dtype=int)
    assert len(s) == nr, (len(s), nr)
    rec[:, I_STATE - 1] = s
if 'poses' in opts:
    p = np.loadtxt(opts['poses'], ndmin=2)
    assert p.shape == (nr, 5), p.shape
    for j, k in enumerate((I_E1, I_E2, I_E3, I_X, I_Y)):
        rec[:, k - 1] = p[:, j]
buf[off:off + nr * nb] = rec.astype('<f4').tobytes()
open(fout, 'wb').write(buf)
print('%s -> %s: %s' % (fin, fout, ' '.join(k for k in ('states', 'poses') if k in opts) or 'copy'))
