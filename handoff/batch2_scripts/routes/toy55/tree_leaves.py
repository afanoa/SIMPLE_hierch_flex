#!/usr/bin/env python3
"""Leaf labels of a data-derived tree for one half: every leaf (a node that did not split, or a child too small to split)
gets an index; writes one label per project row (0 = in no leaf, e.g. junk removed at the root) and the leaf names.
usage: tree_leaves.py HALF_RUN_DIR TAG OUT_LABELS"""
import sys, os, glob, numpy as np
run, tag, out = sys.argv[1:4]; base = '%s/%s' % (run, tag)
nodes = sorted([os.path.basename(d) for d in glob.glob(base + '/r*') if os.path.isfile(d + '/subset.txt')], key=lambda s: [int(x) for x in s.split('.')[1:]])
nrow = len(np.loadtxt(base + '/r/children.txt')) if os.path.isfile(base + '/r/children.txt') else None
leaves = []
for nd in nodes:
    kids = [m for m in nodes if m.rsplit('.', 1)[0] == nd and m != nd]; ch = base + '/' + nd + '/children.txt'
    lab = np.loadtxt(ch).astype(int) if os.path.isfile(ch) else None
    if nrow is None and lab is not None: nrow = len(lab)
    if lab is None or lab.max() <= 1:   # did not split: the node is a leaf
        if not kids: leaves.append((nd, np.loadtxt(base + '/' + nd + '/subset.txt', dtype=int, ndmin=1) - 1))
        continue
    for c in range(1, lab.max() + 1):   # split: children without a directory are leaves
        if '%s.%d' % (nd, c) not in nodes: leaves.append(('%s.%d' % (nd, c), np.where(lab == c)[0]))
L = np.zeros(nrow, int)
for i, (nm, rows) in enumerate(leaves): L[rows] = i + 1
jf = lambda nm: ' junk' if os.path.isfile('%s/%s/junk.flag' % (base, nm)) else ''   # junk leaves of dtree_v11 (root-child junk flag)
np.savetxt(out, L, fmt='%d'); open(out + '.names', 'w').write('\n'.join('%d %s %d%s' % (i + 1, nm, len(r), jf(nm)) for i, (nm, r) in enumerate(leaves)) + '\n')
print('%s: %d leaves: %s' % (run.rstrip('/').split('/')[-1], len(leaves), ' '.join('%s(%d%s)' % (nm, len(r), jf(nm)) for nm, r in leaves)))
