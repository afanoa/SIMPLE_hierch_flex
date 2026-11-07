#!/usr/bin/env python3
"""Per-node view of a running (or finished) hetero_tree on a multi-motion toy: for every node with split/split_children.txt (one child label
per node row, 0 = junk / unassigned) the children's sizes, junk share, majority combination, per-factor state-2 fraction and Cramer's V of
the child label with each factor (genuine rows) -> which motion each split follows, while the tree is still running.
usage: mm_node_score.py TREE_ROOT(1_hetero_tree copy with node_*/rows.txt and node_*/split/split_children.txt) truth=.. factors=.."""
import sys, os, glob, numpy as np
args = dict(a.split('=', 1) if '=' in a else ('tree', a) for a in sys.argv[1:]); T = args['tree']
t = np.loadtxt(args['truth'], ndmin=2); rows = t[:, 0].astype(int); combo = t[:, 1].astype(int); junk = t[:, 2] > 0
F = np.loadtxt(args['factors'], ndmin=2); nf = F.shape[1] - 3; fmap = {int(x[0]): x[1:1 + nf].astype(int) for x in F}
fs = np.zeros((len(rows), nf), int)
for i in np.where(~junk)[0]: fs[i] = fmap[combo[i]]
pos = np.zeros(rows.max() + 1, int); pos[rows] = np.arange(len(rows))
def cramer(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]
    if a.max() == 0 or b.max() == 0: return 0.
    o = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(o, (a, b), 1); e = o.sum(1)[:, None] * o.sum(0)[None] / o.sum()
    return float(np.sqrt(((o - e) ** 2 / np.maximum(e, 1e-12)).sum() / o.sum() / (min(o.shape) - 1)))
ds = sorted(glob.glob(T + '/node_r*/split/split_children.txt'), key=lambda p: (len(p), p))
for p in ds:
    nd = p.split('/node_')[1].split('/')[0]; nr = np.loadtxt(os.path.dirname(os.path.dirname(p)) + '/rows.txt', dtype=int, ndmin=1)
    ch = np.loadtxt(p, dtype=int, ndmin=1)
    if len(ch) != len(nr): print('%s: %d labels for %d rows' % (nd, len(ch), len(nr))); continue
    ix = pos[nr]; g = ~junk[ix] & (ch > 0)
    V = [cramer(ch[g], fs[ix][g, q]) for q in range(nf)]; vj = cramer(ch[ch > 0], junk[ix][ch > 0].astype(int))
    print('%-7s %6d rows  K %d  V %s  junkV %.2f  follows %s' % (nd, len(nr), len(set(ch[ch > 0])), ' '.join('f%d %.2f' % (q + 1, v) for q, v in enumerate(V)), vj,
          ','.join('f%d' % (q + 1) for q, v in enumerate(V) if v >= 0.3) or '-'))
    for c in sorted(set(ch)):
        m = ch == c; mg = m & ~junk[ix]; cb, cn = np.unique(combo[ix][m], return_counts=True); o = np.argsort(-cn)[:3]
        fr = ' '.join('f%d %.2f' % (q + 1, (fs[ix][mg, q] == 2).mean() if mg.any() else np.nan) for q in range(nf))
        print('    child %d: %6d rows junk %.2f  state-2 fraction %s  top %s' % (c, m.sum(), junk[ix][m].mean(), fr, ' '.join('%d(%.2f)' % (cb[i], cn[i] / m.sum()) for i in o)))
