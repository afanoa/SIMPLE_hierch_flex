#!/usr/bin/env python3
"""Whole-tree UMAP coloured by LEAF state (the reproducible cell clusters of the region boards), one panel per family;
a particle that belongs to leaf states of several regions is coloured by the first region in sorted order.
usage: umap_leaves.py RUNDIR OUT.png [DROPPED.txt]  (leaf names f1__1_s11 ... drawn black with red labels)"""
import sys, glob, re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
run, out = sys.argv[1:3]; dropped = set(l.split()[0] for l in open(sys.argv[3]) if l.strip()) if len(sys.argv) > 3 else set(); Z = np.load(run + '/fig_global_embedding.npz', allow_pickle=True); emb, pinds, fam = Z['emb'], Z['pinds'].astype(int), Z['family'].astype(int)
leaf = {}; names = []; dropflag = []
for f in sorted(glob.glob(run + '/fam?/region*/level1/states_clusters.txt')):
    m = re.search(r'fam(\d)/region([^/]+)/', f); st = np.loadtxt(f, dtype=int)
    for s in range(1, st.max() + 1):
        rows = np.where(st == s)[0]; nm = 'f%s %s s%d' % (m.group(1), m.group(2).replace('motion_', 'm').replace('_', 'r'), s); names.append(nm); k = len(names); isdrop = ('f%s_%s_s%02d' % (m.group(1), m.group(2), s)) in dropped; dropflag.append(isdrop)
        leaf.setdefault(int(m.group(1)), []).append((k, nm, set(rows.tolist())))
fams = sorted(leaf); fig, axs = plt.subplots(1, len(fams), figsize=(9 * len(fams), 8.5), squeeze=False)
for ax, F in zip(axs[0], fams):
    mF = fam == F; ax.scatter(emb[mF, 0], emb[mF, 1], s=1, c='lightgrey', alpha=0.4); lab = np.zeros(len(pinds), int); taken = np.zeros(len(pinds), bool)
    for k, nm, rows in leaf[F]:
        m = np.isin(pinds - 1, list(rows)) & mF & ~taken; lab[m] = k; taken |= m
    ks = sorted(set(lab[lab > 0].tolist())); cmap = plt.get_cmap('tab20', 20)
    for i, k in enumerate(ks):
        mk = lab == k; dr = dropflag[k - 1]; ax.scatter(emb[mk, 0], emb[mk, 1], s=6 if dr else 2, color='black' if dr else cmap(i % 20), alpha=0.9 if dr else 0.7, zorder=3 if dr else 2); cx, cy = np.median(emb[mk], 0)
        ax.text(cx, cy, names[k - 1] + (' X' if dr else ''), fontsize=7 if dr else 6, color='red' if dr else 'black', weight='bold' if dr else 'normal', ha='center', va='center', zorder=4, bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.8, lw=0))
    ax.set_title('family %d: %d leaf states; black points / red X labels = the %d removed by the support filter (own half-map FSC > 12 A or < 1000 particles)' % (F, len(ks), sum(dropflag[k - 1] for k in ks)), fontsize=8); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); fig.savefig(out, dpi=120); print('->', out, 'leaf states:', len(names))
