#!/usr/bin/env python3
# pairwise agreement of stage-B labels and split children between production and the replays of one node
# (label match maximised over permutations), plus the truth readouts of each.  usage: pa_agree.py KEY NODE TRUTH
import sys, os, glob, itertools, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent')
import latent_states_v4 as LS
from ppinfo_site import read_dump
A = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor'
key, node, truth = sys.argv[1:4]
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
kind, T = LS.load_truth(truth, node, rows)
runs = {'prod': (os.path.join(A, 'drift', key + '_b.txt'), os.path.join(node, 'split', 'split_children.txt'))}
for r in sorted(glob.glob(os.path.join(A, 'replays', key + '_*'))):
    if os.path.exists(os.path.join(r, 'stageB', 'b_ptcl3D.txt')) and os.path.exists(os.path.join(r, 'split', 'split_children.txt')):
        runs[os.path.basename(r)[len(key) + 1:]] = (os.path.join(r, 'stageB', 'b_ptcl3D.txt'), os.path.join(r, 'split', 'split_children.txt'))
lab, ch = {}, {}
for n, (fb, fc) in runs.items():
    lab[n] = read_dump(fb)['state'].astype(int)
    c = np.loadtxt(fc, dtype=int, ndmin=1)
    ch[n] = c[rows - 1] if len(c) > len(rows) else c

def agree(x, y):
    m = (x > 0) & (y > 0)
    ux, uy = np.unique(x[m]), np.unique(y[m])
    best = 0.
    if len(ux) <= 4 and len(ux) == len(uy):
        for p in itertools.permutations(uy):
            mp = dict(zip(ux, p))
            best = max(best, np.mean(np.array([mp[v] for v in x[m]]) == y[m]))
    return best, LS.nmi(x[m], y[m])

def tline(st):
    if kind == 'closed':
        rd = LS.readout(kind, T, st)
        return 'closed enrich %.3f' % rd['best_enrich']
    rd = LS.readout(kind, T, st)
    s = 'NMI %.3f' % rd['nmi_combo']
    for v, nm in ((18, 'closed'), (101, 'holo')):
        if (T['combo'] == v).sum() > 50:
            pos = T['combo'] == v
            best = 0.
            for u in np.unique(st[st > 0]):
                mk = st == u
                tp = (mk & pos).sum()
                if tp:
                    rec, pur = tp / pos.sum(), tp / mk.sum(); best = max(best, 2 * rec * pur / (rec + pur))
            s += ' %s F1 %.2f' % (nm, best)
    return s

names = list(runs)
print('## %s (%d rows)' % (key, len(rows)))
for n in names:
    k = len(np.unique(ch[n][ch[n] > 0]))
    print('  %-8s stage B %-30s | children %d: %s %s' % (n, tline(lab[n]), k, ' '.join('%d' % (ch[n] == u).sum() for u in np.unique(ch[n][ch[n] > 0])),
          tline(ch[n]) if k > 1 else ''))
print('  pairwise stage-B label agreement (permutation-matched) / NMI:')
for a, b in itertools.combinations(names, 2):
    ag, nm = agree(lab[a], lab[b]); agc, nmc = agree(ch[a], ch[b]) if len(np.unique(ch[a][ch[a] > 0])) == len(np.unique(ch[b][ch[b] > 0])) else (float('nan'), float('nan'))
    print('    %-8s vs %-8s  labels %.3f (NMI %.3f)   children %.3f' % (a, b, ag, nm, agc))
