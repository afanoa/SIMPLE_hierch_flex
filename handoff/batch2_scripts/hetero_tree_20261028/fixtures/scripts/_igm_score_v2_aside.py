#!/usr/bin/env python3
"""IgM toy (armsigm1) groups against the simulation truth: a node's split groups (NODE = node id; split/split_children.txt over the
node's rows) or the final leaves (NODE = leaves; hetero_tree_labels.txt). Per group size and junk share; over the genuine (non-junk)
grouped particles the NMI with the arm-configuration combo and eta^2 of each arm's two angles; per arm angle the 10 % tails (each
motion's two extremes): the best group's tail purity (share of its genuine particles in the tail; 0.10 by chance) and the tail's
recall in that group; an extreme counts as delivered when some group of >= MINPOP holds it at purity >= 0.5. Junk: the share of
the 40,000 junk particles in groups that are >= 50 % junk.
NODE = files=ROWS:LABELS scores a replay directly: ROWS the node's project rows, LABELS one group per project row (0 = none).
usage: igm_score.py TREE_DIR NODE|leaves|files=ROWS:LABELS [TRUTH] [MINPOP]"""
import sys, numpy as np
from pathlib import Path
T, node = Path(sys.argv[1]), sys.argv[2]
truth = Path(sys.argv[3] if len(sys.argv) > 3 else '/mnt/beegfs/elmlund/afan/rnasep_klem/factorial/armsigm1_truth/arms_truth.txt')
minpop = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
cols = open(truth).readline().lstrip('#').split(); t = np.loadtxt(truth, comments='#', ndmin=2)
C = {c: t[:, i] for i, c in enumerate(cols)}; N = int(C['row'].max()); idx = np.zeros(N + 1, int); idx[C['row'].astype(int)] = np.arange(len(t))
lab = np.zeros(N + 1, int)
if node == 'leaves':
    l = np.loadtxt(T / 'hetero_tree_labels.txt', dtype=int, ndmin=1); lab[1:len(l) + 1] = l[:N]
elif node.startswith('files='):
    fr, fl = node[6:].split(':'); rows = np.loadtxt(fr, dtype=int, ndmin=1); l = np.loadtxt(fl, dtype=int, ndmin=1)
    lab[rows] = l[rows - 1] if len(l) >= rows.max() else l
else:
    d = T / ('node_' + node); rows = np.loadtxt(d / 'rows.txt', dtype=int, ndmin=1)
    ch = np.loadtxt(d / 'split' / 'split_children.txt', dtype=int, ndmin=1); ch = ch[rows - 1] if len(ch) != len(rows) else ch
    lab[rows] = ch
r = np.nonzero(lab > 0)[0]; g = lab[r]; ti = idx[r]; junk = C['junk'][ti] > 0; gen = ~junk
def nmi(a, b):
    if len(a) == 0: return float("nan")
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; n = len(a)
    pab = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(pab, (a, b), 1); pab /= n; pa, pb = pab.sum(1), pab.sum(0); nz = pab > 0
    mi = (pab[nz] * np.log(pab[nz] / (pa[:, None] * pb[None, :])[nz])).sum()
    return mi / max(np.sqrt((-(pa[pa > 0] * np.log(pa[pa > 0])).sum()) * (-(pb[pb > 0] * np.log(pb[pb > 0])).sum())), 1e-12)
def eta2(v, gg):
    u, inv = np.unique(gg, return_inverse=True); m = np.bincount(inv, v) / np.bincount(inv); return ((m[inv] - v.mean()) ** 2).sum() / max(((v - v.mean()) ** 2).sum(), 1e-12)
ug, cnt = np.unique(g, return_counts=True); jg = {k: junk[g == k].mean() for k in ug}
print('%s %s: %d grouped particles in %d groups (sizes %s); NMI(groups, combo) %.3f over %d genuine' % (T.parent.name, node, len(r), len(ug),
      ' '.join(str(c) for c in sorted(cnt, reverse=True)[:12]) + (' ..' if len(ug) > 12 else ''), nmi(g[gen], C['combo'][ti][gen]), gen.sum()))
jt = int((C['junk'] > 0).sum()); jin = sum(int((junk & (g == k)).sum()) for k in ug if jg[k] >= 0.5)
print('junk: %d of the %d junk particles in groups that are >= 50 %% junk (%d such groups); junk share per group: %s' % (jin, jt,
      sum(1 for k in ug if jg[k] >= 0.5), ' '.join('%d:%.2f' % (k, jg[k]) for k in ug[:12])))
cb, cbn = np.unique(C['combo'][ti][gen], return_counts=True); best = []
for c, cn in zip(cb, cbn):
    m = C['combo'][ti] == c
    best.append(max(2 * (m & gen & (g == k)).sum() / (cn + (gen & (g == k)).sum()) for k in ug))
print('combos: %d, median size %d; best-group F1 per combo median %.3f, max %.3f (a combo is ~%d particles, groups >= %d: F1 >= 0.5 is out of reach)'
      % (len(cb), int(np.median(cbn)), np.median(best), max(best), int(np.median(cbn)), minpop))
big = [k for k in ug if (gen & (g == k)).sum() >= minpop]; dl = 0; tot = 0
print('%-6s %6s | %-34s | %-34s' % ('angle', 'eta^2', 'low 10 % tail: purity, recall (group)', 'high 10 % tail: purity, recall (group)'))
for a in [c for c in cols if c.startswith('a') and c[1].isdigit()]:
    v = C[a][ti][gen]; gg = g[gen]; lo, hi = np.percentile(C[a][C['junk'] == 0], [10, 90]); out = []
    for tail in (v <= lo, v >= hi):
        tot += 1; ps = [(tail[gg == k].mean(), tail[gg == k].sum() / max(tail.sum(), 1), k) for k in big] or [(0., 0., 0)]
        p, rc, k = max(ps); dl += p >= 0.5; out.append('%.2f, %.2f (g%d)%s' % (p, rc, k, ' DELIVERED' if p >= 0.5 else ''))
    print('%-6s %6.3f | %-34s | %-34s' % (a, eta2(v, gg), out[0], out[1]))
print('extremes delivered (a group >= %d at tail purity >= 0.5): %d of %d' % (minpop, dl, tot))
