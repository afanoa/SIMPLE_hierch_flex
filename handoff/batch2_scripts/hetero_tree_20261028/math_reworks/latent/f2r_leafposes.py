#!/usr/bin/env python3
"""f2r_leafposes.py: every particle's final pose. The root project's ptcl3D dump with the e1 e2 e3 x y tokens of each particle
replaced by those of its leaf node's stage-B dump (node-local rows, mapped through the node's rows.txt); particles without a
leaf (label 0) or without a dump keep the root's poses.
usage: f2r_leafposes.py TREE ROOT_DUMP LEAFDUMP_DIR OUT   (leaf dumps named leaf_<node>_ptcl3D.txt)"""
import os, re, sys
import numpy as np

T, rootdump, leafdir, out = sys.argv[1:5]
lab = np.loadtxt(os.path.join(T, 'hetero_tree_labels.txt'), dtype=int)
leaves = [l.split()[0] for l in open(os.path.join(T, 'hetero_tree_leaves.txt')) if l.strip()]
leaf_node = {j + 1: w.split('_')[0] for j, w in enumerate(leaves)}
lines = open(rootdump).read().splitlines()
assert len(lines) == len(lab), 'root dump %d rows, labels %d' % (len(lines), len(lab))
pat = re.compile(r'\b(e1|e2|e3|x|y)=([-+0-9.eE]+)')
nset = 0; missing = []
for x in sorted(set(leaf_node.values())):
    f = os.path.join(leafdir, 'leaf_%s_ptcl3D.txt' % x)
    rows = np.loadtxt(os.path.join(T, 'node_' + x, 'rows.txt'), dtype=int, ndmin=1)
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        missing.append(x); continue
    ll = open(f).read().splitlines()
    if len(ll) != len(rows):
        print('%s: dump %d rows, node %d; skipped' % (x, len(ll), len(rows))); missing.append(x); continue
    its = [l for l, nx in leaf_node.items() if nx == x]
    want = set((np.where(np.isin(lab, its))[0] + 1).tolist())
    for r, line in zip(rows.tolist(), ll):
        if r in want:
            vals = dict(pat.findall(line))
            lines[r - 1] = pat.sub(lambda m: '%s=%s' % (m.group(1), vals[m.group(1)]), lines[r - 1]); nset += 1
open(out, 'w').write('\n'.join(lines) + '\n')
print('final poses for %d of %d rows (%d assigned to leaves); nodes without a usable dump: %s' % (nset, len(lab), int((lab > 0).sum()), ' '.join(missing) or 'none'))
