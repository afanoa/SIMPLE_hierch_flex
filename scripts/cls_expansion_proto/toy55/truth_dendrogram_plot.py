#!/usr/bin/env python3
"""Dendrogram figure of the reference conformer hierarchy (truth_dendrogram.py output): average linkage on the relative
RMS difference of the 10 A conformer maps; leaves labelled conformer (particles), main clades coloured, bar of populations.
usage: truth_dendrogram_plot.py TRUTH_NPZ OUT_PNG [COLOR_CUT]"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
z = np.load(sys.argv[1]); ids, D, pop = z['ids'], z['D'], z['pop']; cut = float(sys.argv[3]) if len(sys.argv) > 3 else 0.40
Z = linkage(squareform(D, checks=False), 'average')
fig, (ax, axp) = plt.subplots(1, 2, figsize=(11, 13), gridspec_kw={'width_ratios': [4, 1], 'wspace': 0.02})
dn = dendrogram(Z, orientation='left', labels=['conf %d' % k for k in ids], color_threshold=cut, above_threshold_color='#888888', leaf_font_size=9, ax=ax)
ax.axvline(cut, color='#bbbbbb', ls='--', lw=1); ax.set_xlabel('relative RMS difference of the 10 A maps (average linkage)')
ax.set_title('Reference hierarchy of the 55 simulated RNase P conformers', fontsize=12, loc='left')
order = dn['leaves']; y = 5 + 10 * np.arange(len(order)); cols = dn['leaves_color_list']
axp.barh(y, pop[order], height=7, color=cols); axp.set_xlabel('particles (200k set)'); ax.tick_params(axis='y', labelleft=False, labelright=False, length=0)
axp.set_ylim(ax.get_ylim()); axp.yaxis.tick_right(); axp.set_yticks(y); axp.set_yticklabels(['conf %d' % ids[i] for i in order], fontsize=9); axp.tick_params(axis='y', length=0)
for s in ('top', 'right'): ax.spines[s].set_visible(False); axp.spines[s].set_visible(False)
ax.text(cut, ax.get_ylim()[0], ' clades joined below %.2f coloured' % cut, va='bottom', ha='right', fontsize=8, color='#666666')
plt.savefig(sys.argv[2], dpi=110, bbox_inches='tight'); print('png', sys.argv[2], '| clades at the cut:', len(set(cols)) - (1 if '#888888' in cols else 0))
