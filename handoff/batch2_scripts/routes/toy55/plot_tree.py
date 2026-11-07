#!/usr/bin/env python3
"""Dendrogram of a state_merge.py run (merge_tree.txt) with the leaves coloured by their closed-arm occupancy.
usage: plot_tree.py merge_tree.txt arm_states.txt out.png [title]"""
import sys, re, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
mt, armf, out = sys.argv[1:4]; title = sys.argv[4] if len(sys.argv) > 4 else ''
hdr = open(mt).readline(); names = hdr.split('are the leaves ')[1].split(')')[0].split(); N = len(names)
S = np.loadtxt(mt, comments='#'); arm = {}
for l in open(armf):
    if l.strip() and not l.startswith('#'): w = l.split(); arm[str(int(re.findall(r'(\d+)\.mrc$', w[0])[0]))] = float(w[1])
size = {i: 1 for i in range(N)}; Z = []
for st, a, b, r, z, d2, na, nb in S:
    a, b, r = int(a), int(b), int(r); size[r] = size[a] + size[b]; Z.append([a, b, d2, size[r]])
Z = np.array(Z); Z[:, 2] = np.maximum.accumulate(Z[:, 2])   # monotone heights for drawing (complete / average linkage already are)
v = np.array([arm.get(n, np.nan) for n in names]); cm = plt.get_cmap('viridis')
fig, ax = plt.subplots(figsize=(7, max(6, 0.075 * N + 2)))
R = dendrogram(Z, orientation='left', labels=names, leaf_font_size=4 if N > 120 else 7, color_threshold=0, above_threshold_color='0.4', ax=ax)
for t in ax.get_yticklabels():
    x = arm.get(t.get_text(), np.nan); t.set_color(cm((x - 0.1) / 0.85) if np.isfinite(x) else 'k'); t.set_fontweight('bold' if x >= 0.8 else 'normal')
sc = ax.scatter([], [], c=[], cmap=cm, vmin=0.1, vmax=0.95); cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.08); cb.set_label('leaf closed-arm occupancy')
ax.set_xlabel('merge height: relative difference energy of the two branches (cross-half, 40-10 A)'); ax.set_title(title, fontsize=9)
fig.tight_layout(); fig.savefig(out, dpi=170); print('wrote', out)
