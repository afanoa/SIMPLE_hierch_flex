#!/usr/bin/env python3
"""rflex6split node decision from one or more split tests of the node's refined states (splitbench16 output dirs, comma-separated, e.g.
the test at the root project's consensus poses and the test at the node's own refined poses): each test maps every unit (refined state)
to the supported group holding most of its particles; two units share a final group when ANY test puts them together (union-find), so a
split survives only if every test supports it; a unit set aside as junk by any test is a junk leaf; partial units are noted. Per final
group: family (>= MINP particles, depth < MAXD, at least two groups: the node splits) or leaves; one group = the node is one leaf.
usage: rflex6split_decide.py STATES TESTDIR[,TESTDIR2] DEPTH [--minp 50000 --maxd 2]  -> stdout: kind particles members(+) note"""
import argparse, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('states'); ap.add_argument('tests'); ap.add_argument('depth', type=int)
ap.add_argument('--minp', type=int, default=50000); ap.add_argument('--maxd', type=int, default=2); a = ap.parse_args()
lab = np.loadtxt(a.states, dtype=int); us = [int(u) for u in np.unique(lab[lab > 0])]; maps, junk, partial, jm = [], set(), set(), {}
for t in a.tests.split(','):
    ch = np.loadtxt(t + '/children_A.txt', dtype=int) + np.loadtxt(t + '/children_B.txt', dtype=int)
    jf, pf = t + '/junk_units.txt', t + '/partial_units.txt'
    jt = set(int(x) for x in open(jf).read().split()) if os.path.isfile(jf) else set(); junk |= jt
    if os.path.isfile(pf): partial |= set(int(x) for x in open(pf).read().split())
    for l in open(t + '/bench.log', errors='replace'):
        m = re.match(r'\s*junk screen unit (\d+): projection on the screened consensus (\S+) / (\S+), rel. SSNR (\S+) / (\S+)', l)
        if m and int(m.group(1)) not in jm: jm[int(m.group(1))] = 'projection %s/%s rel.SSNR %s/%s (%s)' % (m.groups()[1:] + (os.path.basename(t),))
    to = {}
    for u in us:
        if u in jt: continue
        c = ch[lab == u]; c = c[c > 0]; to[u] = int(np.bincount(c).argmax()) if len(c) else 0
    maps.append(to)
partial -= junk; par = {u: u for u in us if u not in junk}
def find(u):
    while par[u] != u: par[u] = par[par[u]]; u = par[u]
    return u
for to in maps:
    by = {}
    for u, g in to.items():
        if u in par and g > 0: by.setdefault(g, []).append(u)
    for m in by.values():
        for v in m[1:]: par[find(v)] = find(m[0])
comps = {}
for u in par: comps.setdefault(find(u), []).append(u)
groups = sorted((sorted(m) for m in comps.values()), key=lambda m: -int(np.isin(lab, m).sum())); split = len(groups) >= 2
for g, mem in enumerate(groups, 1):
    n = int(np.isin(lab, mem).sum()); pm = [x for x in mem if x in partial]
    print('family' if (split and n >= a.minp and a.depth < a.maxd) else 'leaves', n, '+'.join(map(str, mem)), 'joint group %d of %d (%s)%s' % (g, len(groups), ' / '.join('%d groups' % len(set(v for v in to.values() if v > 0)) for to in maps), '; partial units %s' % '+'.join(map(str, pm)) if pm else ''))
for s in sorted(junk): print('junk', int((lab == s).sum()), s, 'junk screen: ' + jm.get(s, ''))
