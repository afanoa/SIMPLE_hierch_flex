#!/usr/bin/env python3
# the positive-control toy's node readout: per split child (and stage-B state) its size, kind shares (arm continuum, closed,
# holo, junk) and the arm particles' s1 (mean, p10-p90); the s1 AUROC between the two children (1 = a clean cut along the
# continuum) and the best child's closed / holo F1.  usage: pa_qflex.py KEY NODE
import sys, os, glob, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent')
import latent_states_v4  # puts PPINFO (ppinfo_site) on the path
from ppinfo_site import read_dump
A = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor'
TRUTH = '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoyflex/truth_qtoy.txt'
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
key, node = sys.argv[1:3]
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
t = np.loadtxt(TRUTH, ndmin=2)[rows - 1]
C = {c: t[:, COLS.index(c)] for c in COLS}
kind = np.where(C['junk'] > 0, 0, np.where(C['C_closed'] > 0, 2, np.where(C['holo'] > 0, 3, 1)))
s1 = C['s1']

def auroc(a, b):
    x = np.r_[a, b]; r = np.argsort(np.argsort(x)) + 1.
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))

def f1(g, v):
    pos = kind == v; best = 0.
    for u in np.unique(g[g > 0]):
        m = g == u; tp = (m & pos).sum()
        if tp: rec, pur = tp / pos.sum(), tp / m.sum(); best = max(best, 2 * rec * pur / (rec + pur))
    return best

def show(name, g, what):
    us = [u for u in np.unique(g[g > 0])]
    out = '  %-9s %-8s' % (name, what)
    for u in us:
        m = g == u; am = m & (kind == 1)
        out += ' | %d: n %d arm %.2f closed %.2f holo %.2f junk %.2f s1 %.1f [%.1f,%.1f]' % (u, m.sum(), (kind[m] == 1).mean(),
               (kind[m] == 2).mean(), (kind[m] == 3).mean(), (kind[m] == 0).mean(), s1[am].mean(), *np.percentile(s1[am], [10, 90]))
    if len(us) == 2:
        a = s1[(g == us[0]) & (kind == 1)]; b = s1[(g == us[1]) & (kind == 1)]
        au = auroc(a, b); out += ' || s1 AUROC %.3f' % max(au, 1 - au)
    out += ' | closed F1 %.2f holo F1 %.2f' % (f1(g, 2), f1(g, 3))
    print(out)

runs = {'prod': (os.path.join(A, 'drift', key + '_b.txt'), os.path.join(node, 'split', 'split_children.txt'))}
for r in sorted(glob.glob(os.path.join(A, 'replays', key + '_*'))):
    if os.path.exists(os.path.join(r, 'split', 'split_children.txt')):
        runs[os.path.basename(r)[len(key) + 1:]] = (os.path.join(r, 'stageB', 'b_ptcl3D.txt'), os.path.join(r, 'split', 'split_children.txt'))
print('## %s: %d rows; kinds arm %.3f closed %.3f holo %.3f junk %.3f' % (key, len(rows), *[(kind == k).mean() for k in (1, 2, 3, 0)]))
for n, (fb, fc) in runs.items():
    st = read_dump(fb)['state'].astype(int)
    ch = np.loadtxt(fc, dtype=int, ndmin=1); ch = ch[rows - 1] if len(ch) > len(rows) else ch
    show(n, st, 'stageB')
    show(n, ch, 'children')
