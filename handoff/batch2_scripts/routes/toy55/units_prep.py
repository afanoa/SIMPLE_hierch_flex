#!/usr/bin/env python3
"""usage: units_prep.py TREE_DIR OUT_DIR [--labels STATE_LABELS]
Builds a validate_tree3.py layout for the merge-up units of TREE_DIR (hetero_tree_merged.txt = unit per row) in OUT_DIR:
hetero_tree_labels.txt = unit per row (with --labels, each row's state mapped through the tree's leaf -> unit map; 0 stays 0) and
hetero_tree_leaves.txt = one line per unit named r_s<unit> (a flat partition), marked junk when the unit holds a leaf the tree named
junk. Prints each unit's leaves and, for every inner node of the tree, how many units its leaves fall in (1 = collapsed)."""
import os, sys, numpy as np
D, O = sys.argv[1:3]; opt = dict(zip(sys.argv[3::2], sys.argv[4::2]))
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); U = np.loadtxt(D + '/hetero_tree_merged.txt', dtype=int)
LV = [l.split() for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]; names = [w[0] for w in LV]
leaves = [k for k in np.unique(L) if k > 0]; lu = {}
for k in leaves:
    u, c = np.unique(U[L == k], return_counts=True); lu[k] = int(u[np.argmax(c)])
    if len(u) > 1: print('   WARNING leaf %s spans units %s' % (names[k - 1], dict(zip(u.tolist(), c.tolist()))))
units = sorted(set(lu.values()) - {0}); lmap = np.zeros(len(LV) + 1, int)
for k in leaves: lmap[k] = lu[k]
UL = lmap[np.loadtxt(opt['--labels'], dtype=int)] if '--labels' in opt else U
os.makedirs(O, exist_ok=True); np.savetxt(O + '/hetero_tree_labels.txt', UL, fmt='%d')
with open(O + '/hetero_tree_leaves.txt', 'w') as f:
    for u in range(1, max(units) + 1):
        mem = [k for k in leaves if lu[k] == u]
        f.write('r_s%d %s %d %s\n' % (u, 'junk' if any(LV[k - 1][1] == 'junk' for k in mem) else 'leaf', (UL == u).sum(), '+'.join(names[k - 1] for k in mem) or '-'))
print('== %s%s: %d leaves -> %d units' % (D, ' (labels %s)' % opt['--labels'] if '--labels' in opt else '', len(leaves), len(units)))
for u in units: print('   unit %d (%d rows): %s' % (u, (UL == u).sum(), ' + '.join(names[k - 1] for k in leaves if lu[k] == u)))
def path(nm):
    node, _, rest = nm.partition('_'); d = node[1:]; p = ['r' + d[:i] for i in range(len(d) + 1)]
    return p + [nm] if rest.startswith('s') else p
inner = {}
for k in leaves:
    for x in path(names[k - 1])[:-1]: inner.setdefault(x, set()).add(k)
for x in sorted(inner, key=lambda s: (len(s), s)):
    us = sorted({lu[k] for k in inner[x]}); print('   node %-10s %2d leaves -> %d unit%s' % (x, len(inner[x]), len(us), ' COLLAPSED' if len(us) == 1 else 's %s' % us))
