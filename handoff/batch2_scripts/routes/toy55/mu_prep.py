#!/usr/bin/env python3
"""Merge up the tree, one level: the units of every node of the level -> per-half candidate labels for the gold-standard
split test (splitbench8 --pair identity --xem). Level zoom: units = the consolidated zoom sub-states of one region leaf
(names f<F>_<REG>_s<LL>_s<ZZ> in groups.txt, labels per meta row in states.txt). Level region: units = the previous level's
surviving units of one region board (units.txt + units_A/B.txt of mu_collect.py). A node with one unit passes through.
Writes OUT/<node>/cand_{A,B}.txt + labels.txt (cand -> unit, nA, nB) and OUT/nodes.txt (node dir, focus mask, units).
usage: mu_prep.py zoom RUN CONS_DIR MRA MRB OUT  |  mu_prep.py region|family RUN PREV_DIR MRA MRB OUT (env FOCUS_FAM, {F} = family)"""
import sys, os, re, numpy as np
level, run, src, mra, mrb, out = sys.argv[1:7]; os.makedirs(out, exist_ok=True)
MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in (('A', mra), ('B', mrb))}
def focus(board):   # the board's own focus mask (board.txt), the mask the node was split in
    if board.endswith('/family'): return os.environ.get('FOCUS_FAM', 'none').replace('{F}', board.split('/')[0][3:])   # family level: e.g. the family's periphery
    f = run + '/' + board + '/board/board.txt'
    for l in open(f) if os.path.isfile(f) else []:
        if l.startswith('focus_mask'): return l.split()[1]
    return 'none'
if level == 'zoom':
    st = np.loadtxt(src + '/states.txt', dtype=int); U = {}
    for l in open(src + '/groups.txt'):
        if l.startswith('#') or not l.strip(): continue
        s, n, name = l.split()[:3]; m = re.match(r'^f(\d)_(.*)_s(\d+)_s(\d+)$', name)
        if m: U[int(s)] = (int(m.group(1)), m.group(2), int(m.group(3)), name)
    lab = {h: st[MR[h]] for h in 'AB'}; key = lambda u: (U[u][0], U[u][1], U[u][2]); board = lambda k: 'fam%d/zoomregion%s_s%02d' % k
else:   # region level: units of one region board; family level: every unit of one family (across branches)
    U = {}
    for l in open(src + '/units.txt'):
        if l.startswith('#') or not l.strip(): continue
        w = l.split(); U[int(w[0])] = (int(w[1]), w[2], 0, w[6])
    lab = {h: np.loadtxt('%s/units_%s.txt' % (src, h), dtype=int) for h in 'AB'}
    if level == 'family': key = lambda u: (U[u][0],); board = lambda k: 'fam%d/family' % k
    else: key = lambda u: (U[u][0], U[u][1]); board = lambda k: 'fam%d/region%s' % k
nodes = {}
for u in U: nodes.setdefault(key(u), []).append(u)
with open(out + '/nodes.txt', 'w') as fn:
    fn.write('# node_dir focus_mask units (one line per node; a node with one unit passes through)\n')
    for k in sorted(nodes):
        us = sorted(nodes[k]); d = '%s/%s' % (out, board(k).replace('/', '_')); cnt = {u: tuple(int((lab[h] == u).sum()) for h in 'AB') for u in us}
        us = [u for u in us if min(cnt[u]) > 0]
        if not us: continue
        fn.write('%s %s %s\n' % (d, focus(board(k)), ','.join(str(u) for u in us)))
        if len(us) < 2: continue
        os.makedirs(d, exist_ok=True)
        for h in 'AB':
            c = np.zeros(len(lab[h]), int)
            for i, u in enumerate(us): c[lab[h] == u] = i + 1
            np.savetxt('%s/cand_%s.txt' % (d, h), c, fmt='%d')
        with open(d + '/labels.txt', 'w') as f:
            f.write('# cand unit nA nB name\n'); [f.write('%d %d %d %d %s\n' % (i + 1, u, *cnt[u], U[u][3])) for i, u in enumerate(us)]
tested = sum(1 for k in nodes if len(nodes[k]) >= 2)
print('%s level: %d units in %d nodes, %d nodes with >= 2 units to test' % (level, len(U), len(nodes), tested))
