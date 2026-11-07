#!/usr/bin/env python3
# qtoy2 node readout per run: split children (and stage-B states) with size, junk share and the top structure combos;
# purity = the top combo's share of the child's genuine particles; NMI of the children with the combo (genuine particles).
# usage: pa_qtoy.py KEY NODE
import sys, os, glob, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent')
import latent_states_v4 as LS
from ppinfo_site import read_dump
A = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor'
TRUTH = 'qtoy:/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt'
key, node = sys.argv[1:3]
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
kind, T = LS.load_truth(TRUTH, node, rows)
combo, junk = T['combo'].astype(int), T['junk'] > 0
runs = {'prod': (os.path.join(A, 'drift', key + '_b.txt'), os.path.join(node, 'split', 'split_children.txt'))}
for r in sorted(glob.glob(os.path.join(A, 'replays', key + '_*'))):
    if os.path.exists(os.path.join(r, 'split', 'split_children.txt')):
        runs[os.path.basename(r)[len(key) + 1:]] = (os.path.join(r, 'stageB', 'b_ptcl3D.txt'), os.path.join(r, 'split', 'split_children.txt'))

def show(name, g, what):
    us = np.unique(g[g > 0]); out = '  %-8s %-8s %d:' % (name, what, len(us)); pur = []
    for u in us:
        m = (g == u) & ~junk; c, n = np.unique(combo[m], return_counts=True); o = np.argsort(-n)
        pur.append(n[o[0]] / m.sum())
        out += ' | %d: n %d junk %.2f %s' % (u, (g == u).sum(), junk[g == u].mean(), ' '.join('%d(%.2f)' % (c[j], n[j] / m.sum()) for j in o[:2]))
    gm = (g > 0) & ~junk
    out += ' || mean purity %.3f  NMI(genuine) %.3f' % (np.mean(pur), LS.nmi(g[gm], combo[gm]) if len(us) > 1 else 0.)
    print(out)

print('## %s: %d rows, junk %.3f, combos %s' % (key, len(rows), junk.mean(), ' '.join('%d:%d' % (c, n) for c, n in zip(*np.unique(combo[~junk], return_counts=True)) if n > 50)))
for n, (fb, fc) in runs.items():
    st = read_dump(fb)['state'].astype(int)
    ch = np.loadtxt(fc, dtype=int, ndmin=1); ch = ch[rows - 1] if len(ch) > len(rows) else ch
    show(n, st, 'stageB'); show(n, ch, 'children')
