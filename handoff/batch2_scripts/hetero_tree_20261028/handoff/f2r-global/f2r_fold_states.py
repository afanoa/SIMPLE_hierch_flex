#!/usr/bin/env python3
"""f2r_fold_states.py: hetero_tree's stage-A -> stage-B label step on a ptcl3D dump: states of at most MINPOP particles fold into
the largest (fold_small), the labels are compacted to 1..K' (compact; 0 stays 0). Writes one label per row, prints K'.
usage: f2r_fold_states.py stageA_ptcl3D_dump.txt labels_out.txt"""
import os, sys
import numpy as np
sys.path.insert(0, os.environ.get('PPINFO', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo'))
from ppinfo_site import read_dump

MINPOP = 5
st = read_dump(sys.argv[1])['state'].astype(int)
pres = [s for s in np.unique(st) if s > 0]
if len(pres) >= 2:
    big = max(pres, key=lambda s: (st == s).sum())
    for s in pres:
        if (st == s).sum() <= MINPOP:
            st[st == s] = big
pres = [s for s in np.unique(st) if s > 0]
out = np.zeros_like(st)
for j, s in enumerate(pres, 1):
    out[st == s] = j
np.savetxt(sys.argv[2], out, fmt='%d')
print(len(pres))
