#!/usr/bin/env python3
"""Draw a data-derived tree (dtree.sbatch output): one box per node (particles of half A / half B), a stacked bar of its
truth composition (reference clades of the conformer hierarchy cut at --cut, or the conformers themselves when there
are few), parent -> child edges; leaves outlined. usage: dtree_plot.py RUN_A RUN_B TAG GT_A GT_B OUT_PNG [--truth NPZ --cut 0.40]"""
import sys, os, glob, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('runA'); ap.add_argument('runB'); ap.add_argument('tag'); ap.add_argument('gtA'); ap.add_argument('gtB'); ap.add_argument('png')
ap.add_argument('--truth', default=''); ap.add_argument('--cut', type=float, default=0.40); ap.add_argument('--title', default=''); ap.add_argument('--orient', default='auto', help='v: root on top | h: root on the left, leaves stacked (auto: h above 8 leaves)'); a = ap.parse_args()
confA = np.loadtxt(a.gtA, comments='#')[:, 1].astype(int); confB = np.loadtxt(a.gtB, comments='#')[:, 1].astype(int)
nodes = sorted([os.path.basename(d) for d in glob.glob('%s/%s/r*' % (a.runA, a.tag)) if os.path.isfile(d + '/subset.txt')], key=lambda s: [int(x) if x.isdigit() else 0 for x in s.split('.')])
present = sorted(set(confA[confA > 0]) | set(confB[confB > 0]))
if a.truth and len(present) > 6:   # categories = reference clades
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    z = np.load(a.truth); cl = fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=a.cut, criterion='distance'); cmap = dict(zip(z['ids'], cl))
    cats = sorted(set(cmap[k] for k in present)); catof = lambda k: cmap.get(k, -1)
    ref = {24: ('orange', '#ff7f0e'), 5: ('red', '#d62728'), 49: ('green', '#2ca02c'), 1: ('purple', '#9467bd')}   # the reference dendrogram's clade colours, by a member conformer
    named = {cmap[k]: v for k, v in ref.items() if k in cmap}
    catname = {c: '%s clade (%d conf.)' % (named.get(c, ('clade %d' % c, None))[0], sum(1 for k in z['ids'] if cmap[k] == c)) for c in cats}
else:
    cats = present; catname = {c: 'conformer %d' % c for c in cats}; catof = lambda k: k
cols = dict(zip(cats, plt.cm.tab10(np.arange(len(cats)) % 10)))
if a.truth and len(present) > 6: cols.update({c: v[1] for c, v in named.items() if c in cols})
def comp(conf, rows):   # shares of ALL particles: the categories, then junk (grey)
    c = conf[rows]; nj = int((c == 0).sum()); good = c[c > 0]; v = np.array([np.sum([catof(k) == ct for k in good]) for ct in cats], float); return v / max(len(rows), 1), nj / max(len(rows), 1), len(rows)
info = {}; rowsof = {}
for nd in nodes:
    rowsof[nd] = (np.loadtxt('%s/%s/%s/subset.txt' % (a.runA, a.tag, nd), dtype=int, ndmin=1) - 1, np.loadtxt('%s/%s/%s/subset.txt' % (a.runB, a.tag, nd), dtype=int, ndmin=1) - 1)
for nd in list(nodes):   # children too small to split again are leaves without a directory: take them from the node's labels
    fa, fb = '%s/%s/%s/children.txt' % (a.runA, a.tag, nd), '%s/%s/%s/children.txt' % (a.runB, a.tag, nd)
    if not (os.path.isfile(fa) and os.path.isfile(fb)): continue
    la, lb = np.loadtxt(fa).astype(int), np.loadtxt(fb).astype(int); kmax = max(la.max(), lb.max())
    if kmax < 2: continue
    for c in range(1, kmax + 1):
        m = '%s.%d' % (nd, c)
        if m not in nodes: nodes.append(m); rowsof[m] = (np.where(la == c)[0], np.where(lb == c)[0])
nodes = sorted(nodes, key=lambda s: [int(x) if x.isdigit() else 0 for x in s.split('.')])
for nd in nodes: info[nd] = (comp(confA, rowsof[nd][0]), comp(confB, rowsof[nd][1]))
children = {nd: [m for m in nodes if m.rsplit('.', 1)[0] == nd and m != nd] for nd in nodes}
def nleaves(nd): return 1 if not children[nd] else sum(nleaves(c) for c in children[nd])
pos = {}
def place(nd, x0, depth):
    w = nleaves(nd); pos[nd] = (x0 + w / 2., -depth); x = x0
    for c in children[nd]: place(c, x, depth + 1); x += nleaves(c)
place('r', 0., 0); W = nleaves('r'); D = max(-p[1] for p in pos.values()) + 1
horiz = a.orient == 'h' or (a.orient == 'auto' and W > 8)
if horiz:   # depth across, leaves down: x = 2.4 * depth, y = -leaf position
    pos = {nd: (-2.4 * p[1], -p[0]) for nd, p in pos.items()}
    fig, ax = plt.subplots(figsize=(3.2 * D + 2, 0.62 * W + 1.6))
else:
    fig, ax = plt.subplots(figsize=(max(6, 2.2 * W), 1.9 * D + 1.4))
for nd in nodes:
    for c in children[nd]:
        if horiz: ax.plot([pos[nd][0] + 0.45, pos[c][0] - 0.45], [pos[nd][1] - 0.04, pos[c][1] - 0.04], color='#999999', lw=1, zorder=1)
        else: ax.plot([pos[nd][0], pos[c][0]], [pos[nd][1] - 0.18, pos[c][1] + 0.22], color='#999999', lw=1, zorder=1)
for nd in nodes:
    (fa, ja, na), (fb, jb, nb) = info[nd]; x, y = pos[nd]; bw = 0.8; leaf = not children[nd]
    for row, (f, yy) in enumerate(((fa, y + 0.02), (fb, y - 0.10))):
        left = x - bw / 2
        for ct, v in zip(cats, f): ax.barh(yy, v * bw, left=left, height=0.1, color=cols[ct], edgecolor='none', zorder=2); left += v * bw
        jf = ja if row == 0 else jb; ax.barh(yy, jf * bw, left=left, height=0.1, color='#bbbbbb', edgecolor='none', zorder=2)
    ax.add_patch(plt.Rectangle((x - bw / 2 - 0.03, y - 0.17), bw + 0.06, 0.26, fill=False, lw=2.2 if leaf else 0.8, ec='black' if leaf else '#666666', zorder=3))
    ax.text(x, y + 0.13, '%s   A %d / B %d%s' % (nd, na, nb, '' if ja + jb < 0.02 else '  junk %.0f%%' % (50 * (ja + jb))), ha='center', va='bottom', fontsize=8)
if horiz: ax.set_xlim(-0.7, 2.4 * (D - 1) + 0.7); ax.set_ylim(-W + 0.3, 0.6)
else: ax.set_xlim(-0.2, W + 0.2); ax.set_ylim(-D + 0.5, 0.6)
ax.axis('off')
used = [c for i, c in enumerate(cats) if any(info[nd][h][0][i] > 0.005 for nd in nodes for h in (0, 1))]
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=cols[c]) for c in used] + [plt.Rectangle((0, 0), 1, 1, color='#bbbbbb')], labels=[catname[c] for c in used] + ['junk'], loc='lower left', fontsize=8, ncol=min(4, len(cats)), frameon=False, bbox_to_anchor=(0, -0.08))
ax.set_title(a.title or 'data-derived tree (%s): bars = composition of half A (top) and B (bottom); thick boxes = leaves' % a.tag, fontsize=10, loc='left')
plt.tight_layout(); plt.savefig(a.png, dpi=100, bbox_inches='tight'); print('png', a.png, '| nodes', len(nodes), 'leaves', sum(1 for n in nodes if not children[n]))
