#!/usr/bin/env python3
"""usage: mask_vs_truth.py TRUTH_A,TRUTH_B,.. MASK [MASK ..]: Dice and coverage (share of each true site's > 0.5 region inside the mask's > 0.5 region)."""
import sys, os, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
TR = [(os.path.basename(p).replace('.mrc', ''), rd(p) > 0.5) for p in sys.argv[1].split(',')]
for mp in sys.argv[2:]:
    m = rd(mp) > 0.5; print('%s (%d vox): ' % (mp.split('1_hetero_tree/')[-1], m.sum()) + ', '.join('%s Dice %.2f cover %.2f' % (nm, 2. * (m & t).sum() / max(m.sum() + t.sum(), 1), (m & t).sum() / max(t.sum(), 1)) for nm, t in TR))
