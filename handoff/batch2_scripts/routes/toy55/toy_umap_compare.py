#!/usr/bin/env python3
"""Truth-coloured UMAPs of the responsibility vectors of two toy runs on the same particles (80k subsample): is the
conformer structure sharper? Also a silhouette-like separability number per run (15-NN accuracy in UMAP space).
usage: toy_umap_compare.py GT RUN_A:label RUN_B:label OUT.png"""
import sys, numpy as np, umap
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
gt_f, A, B, out = sys.argv[1:5]; gt = {}
for l in open(gt_f):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
runs = [(p.split(':')[0], p.split(':')[1]) for p in (A, B)]; R = [np.load(p + '/klem_resp.npz') for p, _ in runs]; rng = np.random.default_rng(0)
INDEP = len(sys.argv) > 5 and sys.argv[5] == '--independent'
common = np.intersect1d(R[0]['pinds'], R[1]['pinds']) if not INDEP else None; picks = [rng.choice(common, min(80000, len(common)), replace=False)] * 2 if not INDEP else [rng.choice(r['pinds'], min(80000, len(r['pinds'])), replace=False) for r in R]
cmap = plt.get_cmap('tab10'); allc = sorted(set(gt.values()) - {0})
fig, ax = plt.subplots(1, 2, figsize=(18, 8.5))
for j, (r, (p, lab)) in enumerate(zip(R, runs)):
    pick = picks[j]; conf = np.array([gt.get(int(q), 0) for q in pick]); cs = sorted(set(conf[conf > 0]))
    idx = {int(q): i for i, q in enumerate(r['pinds'])}; rows = np.array([idx[int(q)] for q in pick]); X = np.sqrt(r['r'][rows].astype(float))
    U = umap.UMAP(n_neighbors=30, min_dist=0.1, random_state=0).fit_transform(X)
    m = conf > 0; h = np.arange(len(pick)) % 2; acc = []
    for s in (0, 1):
        tr, te = m & (h == s), m & (h != s); nn = cKDTree(U[tr]).query(U[te], k=15)[1]; lab_ = conf[tr][nn]; acc.append(np.mean(np.array([np.bincount(x).argmax() for x in lab_]) == conf[te]))
    ax[j].scatter(U[~m, 0], U[~m, 1], s=1, color='lightgrey', linewidths=0, rasterized=True, label='junk')
    for k, c in enumerate(cs):
        mm = conf == c; ax[j].scatter(U[mm, 0], U[mm, 1], s=1.5, color=cmap(allc.index(c) % 10), alpha=0.5, linewidths=0, rasterized=True, label='conformer %d' % c)
        ax[j].text(np.median(U[mm, 0]), np.median(U[mm, 1]), str(c), fontsize=11, weight='bold', ha='center', bbox=dict(facecolor='white', alpha=0.7, lw=0, pad=1))
    ax[j].set_title('%s: responsibility UMAP, truth colours; 15-NN accuracy in this plane %.2f (chance %.2f)' % (lab, np.mean(acc), np.bincount(conf[m]).max() / m.sum()), fontsize=11); ax[j].set_xticks([]); ax[j].set_yticks([])
ax[0].legend(markerscale=8, fontsize=8, loc='lower left'); plt.tight_layout(); plt.savefig(out, dpi=100); print('figure ->', out)
