#!/usr/bin/env python3
"""usage: set_states_binoris.py IN.simple OUT.simple STATES.txt DUMP.npz
OUT = a copy of IN with ptcl3D's state replaced by STATES (one integer per ptcl3D record), for the projects SIMPLE's selection refuses
(hetero_tree node projects holding only the node's rows: the stk / ptcl2D / ptcl3D counts differ). Project file: a header of 20 segments x
5 int64 (fromto(2), n_bytes_per_record, n_records, first_data_byte; simple_binoris), ptcl3D = segment 6, fixed float32 records. The slots
(I_CORR 3, I_DFX 4, I_E1 8, I_E2 9, I_STATE 23 of simple_defs_ori) are checked against DUMP (dump_ptcl.py of IN) on every record first."""
import shutil, sys, numpy as np
IN, OUT, ST, DP = sys.argv[1:5]
h = np.fromfile(IN, dtype='<i8', count=20 * 5).reshape(20, 5); nb, nr, fdb = (int(x) for x in h[5, 2:])
st = np.loadtxt(ST, dtype=int); D = np.load(DP); assert len(st) == nr == len(D['state']), 'states %d, ptcl3D records %d, dump %d' % (len(st), nr, len(D['state']))
rec = np.memmap(IN, dtype='<f4', mode='r', offset=fdb - 1, shape=(nr, nb // 4))
for k, i in (('corr', 2), ('dfx', 3), ('e1', 7), ('e2', 8), ('state', 22)):
    x, y = rec[:, i].astype(np.float64), D[k]; ok = np.isfinite(y)
    assert ok.sum() > 0.9 * nr and np.allclose(x[ok], y[ok], rtol=1e-4, atol=2e-3), 'slot %d does not hold %s (max |diff| %.3g)' % (i + 1, k, np.abs(x[ok] - y[ok]).max())
del rec; shutil.copyfile(IN, OUT); w = np.memmap(OUT, dtype='<f4', mode='r+', offset=fdb - 1, shape=(nr, nb // 4)); w[:, 22] = st.astype(np.float32); w.flush()
print('%s: %d ptcl3D records (%d bytes each), state set from %s: %s' % (OUT, nr, nb, ST, ', '.join('%d x %d' % (v, (st == v).sum()) for v in np.unique(st))))
