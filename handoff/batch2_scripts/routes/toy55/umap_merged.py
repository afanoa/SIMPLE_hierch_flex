#!/usr/bin/env python3
"""Whole-tree UMAP of a run coloured by the consolidated (merged) state of every particle: left = merged states (grey = not
in any), right = family. usage: umap_merged.py RUNDIR CONSOLIDATED_DIR OUT.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
run, cons, out = sys.argv[1:4]; Z = np.load(run + '/fig_global_embedding.npz', allow_pickle=True); emb, pinds, fam = Z['emb'], Z['pinds'].astype(int), Z['family']
st = np.loadtxt(cons + '/states.txt', dtype=int); lab = st[pinds - 1]; K = lab.max(); groups = [l.split() for l in open(cons + '/groups.txt') if not l.startswith('#')]
fig, axs = plt.subplots(1, 2, figsize=(15, 7)); cmap = plt.get_cmap('tab20', max(K, 1)); o = np.random.default_rng(0).permutation(len(lab))
ax = axs[0]; m = lab == 0; ax.scatter(emb[m, 0], emb[m, 1], s=1, c='lightgrey', alpha=0.4, label='no merged state')
for k in range(1, K + 1):
    mk = lab == k; ax.scatter(emb[mk, 0], emb[mk, 1], s=1.5, color=cmap(k - 1), alpha=0.6); cx, cy = np.median(emb[mk], 0) if mk.any() else (np.nan, np.nan)
    ax.text(cx, cy, str(k), fontsize=11, weight='bold', ha='center', va='center', bbox=dict(boxstyle='round', fc='white', alpha=0.7, lw=0))
ax.set_title('%s: whole-tree UMAP of %d particles coloured by merged state (%d states, numbers = state ids of the merged-state folder)' % (run.rstrip('/').split('/')[-1], len(lab), K), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
ax = axs[1]
for f in sorted(set(fam.tolist())): mf = fam == f; ax.scatter(emb[mf, 0], emb[mf, 1], s=1, alpha=0.5, label='family %d' % f)
ax.legend(markerscale=8, fontsize=9); ax.set_title('the same embedding by family', fontsize=9); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); fig.savefig(out, dpi=120); print('->', out)
for g in groups[:K]: print('state %s: %s particles' % (g[0], g[1]))
