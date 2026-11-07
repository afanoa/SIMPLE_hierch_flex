#!/usr/bin/env python3
"""Score-space embedding of an rflex tree (score_embed.py's embed.npz) against its PLAIN leaf maps, in four panels: the depth-1 family of each
particle's leaf (r1 = the closed family, r5 with its children r51 / r52, the rest), the closed-arm and the open-arm readout of the leaf's plain
map (arm_readout.py on tw_units' raw unit maps, unit = leaf), the old f692r closed particles (states 15+16). Static PNG, light surface.
usage: plot_score_st.py EMBED_DIR ARM_RAW.txt LEAVES.txt OUT.png"""
import sys, re, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
emb, armf, leavesf, png = sys.argv[1:5]; K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'
SURF, INK, INK2, NEUT = '#fcfcfb', '#0b0b0b', '#52514e', '#c9c8c2'; BLUE, ORANGE = '#2a78d6', '#eb6834'
SEQ = LinearSegmentedColormap.from_list('seqblue', ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'])
plt.rcParams.update({'font.size': 9, 'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'font.family': 'DejaVu Sans'})
E = np.load(emb + '/embed.npz', allow_pickle=True); xy, leaf, row, names = E['xy'], E['leaf'], E['row'], [str(x) for x in E['names']]
cnt = [int(l.split()[4]) for l in open(leavesf) if l.strip()]; arm = {}
for l in open(armf):
    if l.startswith('#') or not l.strip(): continue
    p, c, o = l.split()[:3]; arm[int(re.search(r'unit_(\d+)\.mrc', p).group(1))] = (float(c), float(o))
ca = np.array([arm.get(int(l), (np.nan, np.nan))[0] for l in leaf]); oa = np.array([arm.get(int(l), (np.nan, np.nan))[1] for l in leaf])
fam = np.array([1 if names[l - 1].startswith('r1_') else 2 if names[l - 1].startswith('r5') else 3 for l in leaf])
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); oc = ((old == 15) | (old == 16))[row]; kw = dict(s=1.2, linewidths=0, rasterized=True)
def clean(ax): ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
def cbar(sc, ax, lab): cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.01); cb.outline.set_visible(False); cb.set_label(lab, color=INK2)
n1 = sum(c for nm, c in zip(names, cnt) if nm.startswith('r1_')); n5 = sum(c for nm, c in zip(names, cnt) if nm.startswith('r5')); nall = sum(cnt)
fig, axs = plt.subplots(2, 2, figsize=(11, 10.4))
ax = axs[0, 0]; ax.scatter(xy[fam == 3, 0], xy[fam == 3, 1], color=NEUT, alpha=0.35, label='other families (%dk)' % ((nall - n1 - n5) // 1000), **kw)
for k, c, nm in ((2, ORANGE, 'r5 family with r51 / r52 (%dk)' % (n5 // 1000)), (1, BLUE, 'r1, the closed family (%dk)' % (n1 // 1000))):
    m = fam == k; ax.scatter(xy[m, 0], xy[m, 1], color=c, alpha=0.5, label=nm, **kw)
ax.legend(loc='lower left', fontsize=8, frameon=False, markerscale=6); clean(ax); ax.set_title("depth-1 family of each particle's leaf", loc='left', fontsize=10)
ax = axs[0, 1]; o = np.argsort(np.nan_to_num(ca, nan=-1)); sc = ax.scatter(xy[o, 0], xy[o, 1], c=ca[o], cmap=SEQ, vmin=0.3, vmax=1.0, **kw); clean(ax)
cbar(sc, ax, "closed-arm readout of the particle's leaf"); ax.set_title('leaf closed-arm (%d plain leaf maps, max %.2f)' % (len(names), max(v[0] for v in arm.values())), loc='left', fontsize=10)
ax = axs[1, 0]; o = np.argsort(np.nan_to_num(oa, nan=-1)); sc = ax.scatter(xy[o, 0], xy[o, 1], c=oa[o], cmap=SEQ, vmin=0.4, vmax=1.1, **kw); clean(ax)
cbar(sc, ax, "open-arm readout of the particle's leaf"); ax.set_title('leaf open-arm', loc='left', fontsize=10)
ax = axs[1, 1]; ax.scatter(xy[:, 0], xy[:, 1], color=NEUT, alpha=0.3, **kw)
ax.scatter(xy[oc, 0], xy[oc, 1], color=BLUE, s=4, linewidths=0, alpha=0.8, label='old f692r closed states 15+16 (%d of %dk shown)' % (oc.sum(), len(oc) // 1000))
ax.legend(loc='lower left', fontsize=8, frameon=False, markerscale=3); clean(ax); ax.set_title('the old closed particles', loc='left', fontsize=10)
fig.suptitle('%s score space: every particle scored against all %d plain leaf maps (cross-half), PC shares %s; UMAP of %dk'
             % (emb.rstrip('/').split('score_embed_')[-1], len(names), ' '.join('%.2f' % v for v in E['evr'][:5]), len(xy) // 1000), x=0.02, ha='left', fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(png, dpi=160, bbox_inches='tight'); print('wrote ' + png)
hi = np.nan_to_num(ca) >= 0.9; print('UMAP subset %d: old closed %d, of them in leaves with closed-arm >= 0.9: %d (%.1f %%); those leaves hold %d of the subset'
                                     % (len(xy), oc.sum(), (oc & hi).sum(), 100. * (oc & hi).sum() / max(oc.sum(), 1), hi.sum()))
