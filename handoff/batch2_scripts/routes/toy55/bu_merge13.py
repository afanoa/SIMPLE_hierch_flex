#!/usr/bin/env python3
"""bu_merge13.py = bu_merge.py for splitbench13 (junk-aware): units the junk screen set aside (TESTDIR/junk_units.txt, candidate indices)
pass through as their own units and are listed in OUT_PREFIX_junk.txt (new unit ids); listed junk units are never tested again (prep
takes the previous level's junk list). Everything else as bu_merge.py (nested tree cuts G = N/2, N/4, ..., units merged by the majority
of their particles' children; a unit whose particles fall in no child joins the node's largest supported group).
  prep:    bu_merge13.py prep UNITS UNIT_LEAVES JUNK TREE LEVEL_DIR G MRA MRB     (JUNK: file of junk unit ids or '-')
  collect: bu_merge13.py collect UNITS UNIT_LEAVES JUNK LEVEL_DIR OUT_PREFIX MRA MRB
Real data on two half-stacks: MRA / MRB map half rows -> master rows (1-based); candidates and children live in each half's row space."""
import sys, os, numpy as np
def read_ul(f, units):
    if f == '-': return {int(u): [int(u)] for u in np.unique(units) if u > 0}
    return {int(l.split()[0]): [int(x) for x in l.split()[1].split(',')] for l in open(f) if l.strip() and not l.startswith('#')}
def read_junk(f): return set() if f == '-' or not os.path.isfile(f) else {int(x) for x in open(f).read().split()}
def tree_cut(tree, G):
    hdr = open(tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; N = len(names)
    mem = {i: [names[i]] for i in range(N)}; alive = set(range(N)); S = np.atleast_2d(np.loadtxt(tree, comments='#'))
    for st in S[:max(0, N - G)]:
        p, q, r = int(st[1]), int(st[2]), int(st[3]); mem[r] = mem.pop(p) + mem.pop(q); alive -= {p, q}; alive.add(r)
    return [mem[k] for k in sorted(alive)]
mode = sys.argv[1]
if mode == 'prep':
    uf, ulf, jf, tree, ld, G, mra, mrb = sys.argv[2:10]; G = int(G); os.makedirs(ld, exist_ok=True)
    units = np.loadtxt(uf, dtype=int); ul = read_ul(ulf, units); l2u = {lf: u for u, ls in ul.items() for lf in ls}; JK = read_junk(jf)
    MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in zip('AB', (mra, mrb))}; ntest = 0
    with open(ld + '/nodes.txt', 'w') as fn:
        fn.write('# node_dir units (a node with one unit passes through; junk units pass through)\n')
        for gi, g in enumerate(tree_cut(tree, G), 1):
            us = sorted(set(l2u[lf] for lf in g if lf in l2u)); us = [u for u in us if (units[MR['A']] == u).any() and (units[MR['B']] == u).any()]
            for k, u in enumerate([u for u in us if u in JK]): fn.write('%s/n%03d_j%d %d\n' % (ld, gi, k, u))
            us = [u for u in us if u not in JK]
            if not us: continue
            d = '%s/n%03d' % (ld, gi); fn.write('%s %s\n' % (d, ','.join(map(str, us))))
            if len(us) < 2: continue
            os.makedirs(d, exist_ok=True); ntest += 1
            for h in 'AB':
                c = np.zeros(len(MR[h]), int)
                for i, u in enumerate(us): c[units[MR[h]] == u] = i + 1
                np.savetxt('%s/cand_%s.txt' % (d, h), c, fmt='%d')
            with open(d + '/labels.txt', 'w') as f: f.write('# cand unit nA nB\n'); [f.write('%d %d %d %d\n' % (i + 1, u, (units[MR['A']] == u).sum(), (units[MR['B']] == u).sum())) for i, u in enumerate(us)]
    print('cut into %d groups: %d nodes to test (%d units, %d junk passed through)' % (G, ntest, len(ul), len(JK)))
else:
    uf, ulf, jf, ld, outp, mra, mrb = sys.argv[2:9]; units = np.loadtxt(uf, dtype=int); ul = read_ul(ulf, units); JK = read_junk(jf)
    MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in zip('AB', (mra, mrb))}; new = np.zeros_like(units); nl = {}; nid = 0; merged = 0; newjunk = []
    for l in open(ld + '/nodes.txt'):
        if l.startswith('#') or not l.strip(): continue
        d, us = l.split(); us = [int(u) for u in us.split(',')]
        ok = len(us) >= 2 and all(os.path.isfile('%s/children_%s.txt' % (d, h)) for h in 'AB')
        if not ok:
            if len(us) >= 2: print('%s: test missing, %d units passed through' % (os.path.basename(d), len(us)))
            for u in us:
                nid += 1; new[units == u] = nid; nl[nid] = list(ul[u])
                if u in JK: newjunk.append(nid)
            continue
        jc = [int(x) for x in open(d + '/junk_units.txt').read().split()] if os.path.isfile(d + '/junk_units.txt') else []
        ju = [us[c - 1] for c in jc if 1 <= c <= len(us)]
        for u in ju: nid += 1; new[units == u] = nid; nl[nid] = list(ul[u]); newjunk.append(nid)
        ch = np.zeros(len(units), int); ch[MR['A']] = np.loadtxt(d + '/children_A.txt', dtype=int); ch[MR['B']] = np.loadtxt(d + '/children_B.txt', dtype=int)
        rest = [u for u in us if u not in ju]; to = {}
        for u in rest:
            c = ch[units == u]; c = c[c > 0]; to[u] = int(np.bincount(c).argmax()) if len(c) else 0
        big = max(set(to.values()) - {0}, key=lambda g: int((ch == g).sum()), default=1); to = {u: (g if g > 0 else big) for u, g in to.items()}
        for g in sorted(set(to.values())):
            mem = [u for u in rest if to[u] == g]; nid += 1; nl[nid] = sorted(x for u in mem for x in ul[u])
            for u in mem: new[units == u] = nid
            merged += len(mem) - 1
    np.savetxt(outp + '_units.txt', new, fmt='%d'); open(outp + '_junk.txt', 'w').write(''.join('%d\n' % j for j in newjunk))
    with open(outp + '_leaves.txt', 'w') as f: f.write('# unit leaves\n'); [f.write('%d %s\n' % (u, ','.join(map(str, ls)))) for u, ls in nl.items()]
    print('%d units -> %d units (%d merges, %d junk units)' % (len(ul), nid, merged, len(newjunk)))
