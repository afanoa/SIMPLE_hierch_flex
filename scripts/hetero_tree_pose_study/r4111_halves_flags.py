#!/usr/bin/env python3
"""r4111_halves_flags.py NODE_DIR WORK_DIR: from the node's stage-B dump (WORK_DIR/b_dump.txt, node-local rows) and its rows.txt
(node-local -> project row), split the leaf r4111_s1's particles at their median stage-B corr and write two 0/1 flag files over
the master project's rows (WORK_DIR/nr.txt rows): selA.txt = top half, selB.txt = bottom half."""
import sys
import numpy as np
nd, w = sys.argv[1], sys.argv[2]
nr = int(open(w + '/nr.txt').read())
rows = np.loadtxt(nd + '/rows.txt', dtype=int)
leaf = set(np.loadtxt(nd + '/leaf_r4111_s1_rows.txt', dtype=int).tolist())
corr, st = [], []
for ln in open(w + '/b_dump.txt'):
    d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
    corr.append(float(d['corr'])); st.append(int(float(d['state'])))
corr, st = np.array(corr), np.array(st)
assert len(corr) == len(rows), (len(corr), len(rows))
inleaf = np.array([r in leaf for r in rows]) & (st > 0)
med = np.median(corr[inleaf])
top = inleaf & (corr >= med); bot = inleaf & (corr < med)
fa = np.zeros(nr + 1, int); fb = np.zeros(nr + 1, int)
fa[rows[top]] = 1; fb[rows[bot]] = 1
np.savetxt(w + '/selA.txt', fa[1:], fmt='%d'); np.savetxt(w + '/selB.txt', fb[1:], fmt='%d')
print('leaf rows in node %d of %d; median corr %.4f; top %d (corr p50 %.4f) bottom %d (corr p50 %.4f); flags over %d master rows' % (
    inleaf.sum(), len(rows), med, top.sum(), np.median(corr[top]), bot.sum(), np.median(corr[bot]), nr))
