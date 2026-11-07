#!/usr/bin/env python3
"""Merge decisions along the bottom-up tree (generic merge-up, any data set). Units start as the leaves (per-row labels); the leaf tree
(state_merge.py merge_tree.txt) is cut into nested groups per level (G = N/2, N/4, ..., 1); in every group holding >= 2 units, the units
are tested together by splitbench8 on independent halves (the rows' even/odd flag, identity pairing, no re-assignment, no blend pruning):
the groups the halves support as distinct survive (K*), the rest merge. The survivors are the next level's units.
  prep:    bu_merge.py prep UNITS UNIT_LEAVES EO TREE LEVEL_DIR G   -> LEVEL_DIR/nodes.txt, LEVEL_DIR/n<i>/{cand_A,cand_B,labels}.txt
  collect: bu_merge.py collect UNITS UNIT_LEAVES LEVEL_DIR OUT_PREFIX -> OUT_PREFIX_units.txt (per row), OUT_PREFIX_leaves.txt
UNITS: per-row unit labels (0 = none); UNIT_LEAVES: 'unit leaf,leaf,...' per line ('-' at level 0 = every unit is its own leaf).
Real data on two half-stacks: append MRA MRB (row maps half row -> master row, 1-based) to prep and collect; the candidate and
children files are then in each half's row space and EO is ignored (the halves are the stacks)."""
import sys, os, numpy as np
def read_ul(f, units):
    if f == '-': return {int(u): [int(u)] for u in np.unique(units) if u > 0}
    return {int(l.split()[0]): [int(x) for x in l.split()[1].split(',')] for l in open(f) if l.strip() and not l.startswith('#')}
def tree_cut(tree, G):   # leaf groups of the tree after N - G merges (nested for decreasing G)
    hdr = open(tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; N = len(names)
    mem = {i: [names[i]] for i in range(N)}; alive = set(range(N))
    S = np.atleast_2d(np.loadtxt(tree, comments='#'))
    for st in S[:max(0, N - G)]:
        p, q, r = int(st[1]), int(st[2]), int(st[3]); mem[r] = mem.pop(p) + mem.pop(q); alive -= {p, q}; alive.add(r)
    return [mem[k] for k in sorted(alive)]
mode = sys.argv[1]
if mode == 'prep':
    uf, ulf, eof, tree, ld, G = sys.argv[2:8]; G = int(G); os.makedirs(ld, exist_ok=True); MR = None
    units = np.loadtxt(uf, dtype=int); ul = read_ul(ulf, units); l2u = {lf: u for u, ls in ul.items() for lf in ls}
    if len(sys.argv) > 9: MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in zip('AB', sys.argv[8:10])}; eo = np.full(len(units), -1); eo[MR['A']] = 0; eo[MR['B']] = 1
    else: eo = np.loadtxt(eof, dtype=int)
    groups = tree_cut(tree, G); ntest = 0
    with open(ld + '/nodes.txt', 'w') as fn:
        fn.write('# node_dir units (a node with one unit passes through)\n')
        for gi, g in enumerate(groups, 1):
            us = sorted(set(l2u[lf] for lf in g if lf in l2u)); us = [u for u in us if ((units == u) & (eo == 0)).sum() > 0 and ((units == u) & (eo == 1)).sum() > 0]
            if not us: continue
            d = '%s/n%03d' % (ld, gi); fn.write('%s %s\n' % (d, ','.join(map(str, us))))
            if len(us) < 2: continue
            os.makedirs(d, exist_ok=True); ntest += 1
            for h, e in (('A', 0), ('B', 1)):
                if MR is not None:   # half-stack row space
                    c = np.zeros(len(MR[h]), int)
                    for i, u in enumerate(us): c[units[MR[h]] == u] = i + 1
                else:
                    c = np.zeros(len(units), int)
                    for i, u in enumerate(us): c[(units == u) & (eo == e)] = i + 1
                np.savetxt('%s/cand_%s.txt' % (d, h), c, fmt='%d')
            with open(d + '/labels.txt', 'w') as f: f.write('# cand unit nA nB\n'); [f.write('%d %d %d %d\n' % (i + 1, u, ((units == u) & (eo == 0)).sum(), ((units == u) & (eo == 1)).sum())) for i, u in enumerate(us)]
    print('cut into %d groups: %d nodes to test (%d units)' % (len(groups), ntest, len(ul)))
else:
    uf, ulf, ld, outp = sys.argv[2:6]; units = np.loadtxt(uf, dtype=int); ul = read_ul(ulf, units); new = np.zeros_like(units); nl = {}; nid = 0; merged = 0
    MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in zip('AB', sys.argv[6:8])} if len(sys.argv) > 7 else None
    for l in open(ld + '/nodes.txt'):
        if l.startswith('#') or not l.strip(): continue
        d, us = l.split(); us = [int(u) for u in us.split(',')]
        ok = len(us) >= 2 and all(os.path.isfile('%s/children_%s.txt' % (d, h)) for h in 'AB')
        if not ok:
            if len(us) >= 2: print('%s: test missing, %d units passed through' % (os.path.basename(d), len(us)))
            for u in us: nid += 1; new[units == u] = nid; nl[nid] = list(ul[u])
            continue
        if MR is not None:   # children per half-stack row -> master rows
            ch = np.zeros(len(units), int); ch[MR['A']] = np.loadtxt(d + '/children_A.txt', dtype=int); ch[MR['B']] = np.loadtxt(d + '/children_B.txt', dtype=int)
        else: ch = np.loadtxt(d + '/children_A.txt', dtype=int) + np.loadtxt(d + '/children_B.txt', dtype=int)   # disjoint rows (even / odd)
        to = {}
        for u in us:
            c = ch[units == u]; c = c[c > 0]; to[u] = int(np.bincount(c).argmax()) if len(c) else 0
        for g in sorted(set(to.values())):
            mem = [u for u in us if to[u] == g]; nid += 1; nl[nid] = sorted(x for u in mem for x in ul[u])
            for u in mem: new[units == u] = nid
            merged += len(mem) - 1
    np.savetxt(outp + '_units.txt', new, fmt='%d')
    with open(outp + '_leaves.txt', 'w') as f: f.write('# unit leaves\n'); [f.write('%d %s\n' % (u, ','.join(map(str, ls)))) for u, ls in nl.items()]
    print('%d units -> %d units (%d merges)' % (len(ul), nid, merged))
