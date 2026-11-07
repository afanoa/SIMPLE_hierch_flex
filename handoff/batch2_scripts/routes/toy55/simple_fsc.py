#!/usr/bin/env python3
"""usage: simple_fsc.py SERIES_DIR ..
SIMPLE's own even / odd FSC per bin of a motion series: SERIES_DIR/work/fsc_stateKK.bin (float32: shell count, then one value per
shell 1..n/2), the box and pixel size from SERIES_DIR/bin1.mrc; prints the FSC 0.5 / 0.143 resolution (first shell below)."""
import sys, numpy as np
for D in sys.argv[1:]:
    h = np.frombuffer(open(D + '/bin1.mrc', 'rb').read(1024), np.int32); hf = h.view(np.float32); n = int(h[0]); smpd = float(hf[10]) / n; out = []
    for k in range(1, 6):
        a = np.fromfile('%s/work/fsc_state%02d.bin' % (D, k), np.float32); c = a[1:1 + int(a[0])]
        r = [n * smpd / next((i + 1 for i in range(len(c)) if c[i] < t), len(c)) for t in (0.5, 0.143)]; out.append('bin%d %.2f / %.2f' % (k, r[0], r[1]))
    print('simple FSC 0.5 / 0.143 (A): %s' % ' | '.join(out))
