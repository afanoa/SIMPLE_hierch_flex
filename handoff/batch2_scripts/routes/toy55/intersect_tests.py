#!/usr/bin/env python3
"""Join of two split tests of the same units (e.g. consensus poses and the node's refined poses): each test maps every unit to the
supported group holding most of its particles; two units end in the same final group when EITHER test puts them together (union-find),
so a split survives only if both tests support it; a unit is junk when either test set it aside. Prints the final groups (units,
particles, and with --gt their truth composition and NMI(conformer; group)).
usage: intersect_tests.py UNITS TESTDIR1 TESTDIR2 [--gt gt_particles.txt] [--out groups.txt]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('units'); ap.add_argument('t1'); ap.add_argument('t2'); ap.add_argument('--gt', default=''); ap.add_argument('--out', default=''); a = ap.parse_args()
lab = np.loadtxt(a.units, dtype=int); us = [int(u) for u in np.unique(lab[lab > 0])]
def rd(t):
    ch = np.loadtxt(t + '/children_A.txt', dtype=int) + np.loadtxt(t + '/children_B.txt', dtype=int)
    jf = t + '/junk_units.txt'; junk = set(int(x) for x in open(jf).read().split()) if os.path.isfile(jf) else set()
    pf = t + '/partial_units.txt'; part = set(int(x) for x in open(pf).read().split()) if os.path.isfile(pf) else set()
    to = {}
    for u in us:
        if u in junk: continue
        c = ch[lab == u]; c = c[c > 0]; to[u] = int(np.bincount(c).argmax()) if len(c) else 0
    return to, junk, part
(t1, j1, p1), (t2, j2, p2) = rd(a.t1), rd(a.t2); junk = j1 | j2; part = (p1 | p2) - junk
par = {u: u for u in us if u not in junk}
def find(u):
    while par[u] != u: par[u] = par[par[u]]; u = par[u]
    return u
for to in (t1, t2):
    by = {}
    for u, g in to.items():
        if u in par and g > 0: by.setdefault(g, []).append(u)
    for g, m in by.items():
        for v in m[1:]: par[find(v)] = find(m[0])
comps = {}
for u in par: comps.setdefault(find(u), []).append(u)
groups = sorted(comps.values(), key=lambda m: -int(np.isin(lab, m).sum()))
conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int)[:len(lab)] if a.gt else None
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); return I / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
fl = np.zeros(len(lab), int); lines = []
for i, m in enumerate(groups, 1):
    rows = np.isin(lab, m); fl[rows] = i; s = 'group %d: units %s, %d particles%s' % (i, '+'.join(map(str, sorted(m))), rows.sum(), ' (partial units %s)' % '+'.join(map(str, sorted(set(m) & part))) if set(m) & part else '')
    if conf is not None:
        c = conf[rows]; cc = np.bincount(c[c > 0]) if (c > 0).any() else np.zeros(1, int); top = np.argsort(-cc)[:3]; sw = c[(c >= 200) & (c <= 230)] - 200
        s += ' | junk %.2f, %s%s' % ((c == 0).mean(), ' '.join('%d:%.2f' % (t, cc[t] / len(c)) for t in top if cc[t] > 0), ', swing %.1f+-%.1f A (%.0f%%)' % (sw.mean(), sw.std(), 100 * len(sw) / len(c)) if len(sw) > 50 else '')
    lines.append(s)
for i, u in enumerate(sorted(junk), len(groups) + 1):
    rows = lab == u; fl[rows] = i; lines.append('junk unit %d: %d particles%s' % (u, rows.sum(), ' | junk %.2f' % (conf[rows] == 0).mean() if conf is not None else ''))
print('joint decision: %d groups + %d junk units (test 1: %d groups, test 2: %d groups)' % (len(groups), len(junk), len(set(t1.values())), len(set(t2.values()))))
for l in lines: print('  ' + l)
if conf is not None: on = fl > 0; print('  NMI(conformer; unit) %.3f -> joint groups + junk units %.3f' % (nmi(conf[lab > 0], lab[lab > 0]), nmi(conf[on], fl[on])))
if a.out: open(a.out, 'w').write(''.join('%s\n' % '+'.join(map(str, sorted(m))) for m in groups))
