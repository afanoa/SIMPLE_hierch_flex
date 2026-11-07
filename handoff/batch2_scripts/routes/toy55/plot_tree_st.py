#!/usr/bin/env python3
"""Construction tree of an rflex run as a horizontal dendrogram: root -> depth-1 families -> deeper nodes -> leaves (aligned right, DFS order,
children sorted by their subtree's most closed leaf), leaves coloured by the closed-arm readout of their delivered map, with a table per leaf:
file (unit_NN.mrc = leaf NN of leaves.txt), leaf, particles, closed-arm, open-arm, share of the old f692r closed states 15+16, a readout note.
Dark lines = paths to leaves with closed-arm >= 0.9. Static PNG, light surface.  usage: plot_tree_st.py RFLEX_DIR ARM.txt OUT.png"""
import sys, re, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
R, armf, png = sys.argv[1:4]; K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; RT = R.rstrip('/').split('rflex_')[-1]
SURF, INK, INK2, GRID, NEUT = '#fcfcfb', '#0b0b0b', '#52514e', '#d6d5d0', '#c9c8c2'
SEQ = LinearSegmentedColormap.from_list('seqblue', ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'])
plt.rcParams.update({'font.size': 8, 'text.color': INK, 'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'font.family': 'DejaVu Sans'})
L = [l.split() for l in open(R + '/leaves.txt') if l.strip()]; names = [l[0] for l in L]; NL = len(L); idx = {nm: i for i, nm in enumerate(names)}; arm = {}
for l in open(armf):
    if l.startswith('#') or not l.strip(): continue
    p, c, o = l.split()[:3]; arm[int(re.search(r'unit_(\d+)\.mrc', p).group(1))] = (float(c), float(o))
lab = np.loadtxt(R + '/final_labels.txt', dtype=int); old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); assert len(lab) == len(old)
oc = (old == 15) | (old == 16); n = np.bincount(lab, minlength=NL + 1)[1:]; nc = np.bincount(lab[oc], minlength=NL + 1)[1:]; TC = int(nc.sum())
Kin, Kst = {}, {}   # per node: states refined (K) and states the split test kept (K*), last occurrence in the driver logs
for f in sorted(glob.glob(K_ + '/klem_rflex1*.out')):
    for l in open(f, errors='replace'):
        if 'rflex %s]' % RT not in l: continue
        m = re.search(r'node (r\d*) \(depth \d+, \d+ particles\): flex_pca on \d+, (\d+) states', l) or re.search(r'node (r\d*): split test K\* (\d+)', l)
        if m: (Kin if 'flex_pca' in l else Kst)[m.group(1)] = int(m.group(2))
children = {}
for nm in names:   # leaf r51_s5 -> node r51 -> r5 -> r (a node's parent drops its last digit)
    u = nm.split('_')[0]; children.setdefault(u, []).append(nm)
    while u != 'r':
        p = u[:-1]; children.setdefault(p, [])
        if u not in children[p]: children[p].append(u)
        u = p
def stats(u):   # particles, old closed, most closed leaf of the subtree
    if u in idx: i = idx[u]; return n[i], nc[i], arm[i + 1][0]
    s = [stats(c) for c in children[u]]; return sum(x[0] for x in s), sum(x[1] for x in s), max(x[2] for x in s)
for u in children: children[u].sort(key=lambda c: -stats(c)[2])
xs, ys, cur = {}, {}, [0.]
def place(u, d):
    if u in idx: xs[u] = 3.; ys[u] = cur[0]; cur[0] -= 1.; return
    for c in children[u]: place(c, d + 1)
    xs[u] = float(d); ys[u] = np.mean([ys[c] for c in children[u]])
place('r', 0)
hot = {nm for nm in names if arm[idx[nm] + 1][0] >= 0.9}; onpath = set(hot)
for nm in hot:
    u = nm.split('_')[0]
    while True:
        onpath.add(u)
        if u == 'r': break
        u = u[:-1]
fig, ax = plt.subplots(figsize=(13, 2.2 + 0.2 * NL))
for u, ch in children.items():   # elbows: vertical at the parent's x, horizontal to each child
    yc = [ys[c] for c in ch]; ax.plot([xs[u], xs[u]], [min(yc), max(yc)], color=GRID, lw=1.0, zorder=1)
    for c in ch:
        on = c in onpath; ax.plot([xs[u], xs[c]], [ys[c], ys[c]], color=INK if on else GRID, lw=2.0 if on else 1.0, zorder=2 if on else 1)
        if on: ax.plot([xs[u], xs[u]], [ys[u], ys[c]], color=INK, lw=2.0, zorder=2)
NT = stats('r')[0]
for u in children:
    N_, C_, _ = stats(u); ax.scatter(xs[u], ys[u], s=30 + 250 * np.sqrt(N_ / NT), color=NEUT, edgecolor=SURF, linewidth=1.5, zorder=3)
    kk = '\n%d -> %d states kept' % (Kin.get(u, min(8, max(2, N_ // 5000))), Kst[u]) if u in Kst else ''
    ax.text(xs[u] - 0.06, ys[u] + 0.35, '%s  %.1fk%s\n%.0f %% of old closed' % (u, N_ / 1000., kk, 100. * C_ / max(TC, 1)), ha='right', va='bottom', fontsize=7.5,
            color=INK if u in onpath else INK2, fontweight='bold' if u in onpath else 'normal', zorder=4)
cols = [(3.22, 'file', 'left'), (3.72, 'leaf', 'left'), (4.6, 'particles', 'right'), (5.1, 'closed\narm', 'right'), (5.55, 'open\narm', 'right'),
        (6.1, 'old\nclosed', 'right'), (6.25, 'readout note', 'left')]
for x, t, ha in cols: ax.text(x, 1.0, t, ha=ha, va='bottom', fontsize=7.5, color=INK2, fontweight='bold')
nmax = n.max()
for nm in names:
    i = idx[nm]; c, o = arm[i + 1]; y = ys[nm]
    ax.scatter(3., y, s=18 + 110 * n[i] / nmax, color=SEQ(min(c, 1.)), edgecolor=SURF, linewidth=1.0, zorder=3)
    note = 'closed' if c >= 0.9 and o <= 0.6 else 'closed and open regions both filled' if c >= 0.9 else 'arm in neither deposited position' if c < 0.25 and o < 0.65 else ''
    kw = dict(va='center', fontsize=7.5, color=INK if c >= 0.6 else INK2, fontweight='bold' if c >= 0.9 else 'normal', zorder=4)
    for (x, _, ha), t in zip(cols, ('unit_%02d' % (i + 1), nm, '{:,}'.format(int(n[i])), '%.2f' % c, '%.2f' % o, '%.1f %%' % (100. * nc[i] / max(n[i], 1)), note)):
        ax.text(x, y, t, ha=ha, **kw)
ax.set_xlim(-0.75, 8.1); ax.set_ylim(cur[0] - 0.2, 2.9); ax.axis('off')
cax = ax.inset_axes([0.0, -0.012, 0.3, 0.007]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, 1)), cax=cax, orientation='horizontal')
cb.outline.set_visible(False); cb.set_label('leaf closed-arm (relative to the deposited closed conformer)', color=INK2, fontsize=7.5); cb.ax.tick_params(labelsize=7, colors=INK2)
ax.set_title('%s construction tree: %s particles -> %d families -> %d leaves. Each node refines K = min(8, n / 5000) states, the split test merges the ones the halves do not\n'
             'reproduce ("K -> kept"); a kept group becomes a node only with >= 15,000 particles (else a leaf). '
             'Leaf readouts: delivered plain maps (unit_NN.mrc).\nOld closed = share of the old f692r closed states 15+16 (%s particles). '
             'Dark = paths to closed-arm >= 0.9' % (RT, '{:,}'.format(int(NT)), len(children['r']), NL, '{:,}'.format(TC)), loc='left', fontsize=10, color=INK)
fig.savefig(png, dpi=160, bbox_inches='tight'); print('wrote %s (%d leaves, %d nodes, old closed %d)' % (png, NL, len(children), TC))
