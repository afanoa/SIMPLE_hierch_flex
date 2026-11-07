#!/usr/bin/env python3
"""Per-particle poses at the end of a refined tree (dtree_v13): every leaf's particles take the poses of the deepest
refinement that covered them: the leaf node's own refinement (NODE/refined/poses_full.txt), else what the node
inherited (NODE/poses_full.txt), else its parent's; rows in no leaf keep BASE. usage: merge_leaf_poses.py RUN_HALF TAG LEAVES BASE OUT"""
import sys, os, numpy as np
run, tag, leavesf, base, out = sys.argv[1:6]; T = os.path.join(run, tag)
L = np.loadtxt(leavesf).astype(int); names = {int(l.split()[0]): l.split()[1] for l in open(leavesf + '.names') if l.strip()}
lines = open(base).read().splitlines(); src_of = {}
def source(nm):
    while nm:
        for f in ('refined/poses_full.txt', 'poses_full.txt'):
            p = os.path.join(T, nm, f)
            if os.path.exists(p): return os.path.realpath(p)
        nm = nm.rsplit('.', 1)[0] if '.' in nm else ''
    return base
for k, nm in names.items(): src_of.setdefault(source(nm), []).append(k)
for src, ks in src_of.items():
    rows = np.where(np.isin(L, ks))[0]
    if src == os.path.realpath(base) or not len(rows): continue
    S = open(src).read().splitlines(); assert len(S) >= len(lines), (src, len(S), len(lines))
    for i in rows: lines[i] = S[i]
    print('%d rows from %s' % (len(rows), src.replace(run, '')), flush=True)
open(out, 'w').write('\n'.join(lines) + '\n'); print('merged poses -> %s (%d rows)' % (out, len(lines)))
