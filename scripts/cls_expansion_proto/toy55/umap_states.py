#!/usr/bin/env python3
"""Per-particle UMAP from the cross-half state memberships saved by state_graph.py (soft posteriors of every particle over
the final states, scored against the OTHER half's maps, Hellinger-transformed): left coloured by the assigned state, right
by the truth (top conformers, junk grey) when known. Islands here come only from the data's own state evidence; particles
whose evidence is split between states fall on the bridges. usage: umap_states.py GRAPH_NPZ OUT.png [--nsub 60000]"""
import sys, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import umap
ap = argparse.ArgumentParser(); ap.add_argument('npz'); ap.add_argument('png'); ap.add_argument('--nsub', type=int, default=60000); ap.add_argument('--title', default=''); a = ap.parse_args()
Z = np.load(a.npz); W = np.r_[Z['WA'], Z['WB']].astype(np.float32); own = np.r_[Z['ownA'], Z['ownB']]; conf = np.r_[Z['confA'], Z['confB']]; ids = list(Z['ids'])
ok = own >= 0; W, own, conf = W[ok], own[ok], conf[ok]; rng = np.random.default_rng(0); sel = rng.choice(len(W), min(a.nsub, len(W)), replace=False); W, own, conf = W[sel], own[sel], conf[sel]
E = umap.UMAP(n_neighbors=30, min_dist=0.3, metric='euclidean', random_state=0).fit_transform(np.sqrt(np.maximum(W, 0)))
asg = W.argmax(1); pal = plt.cm.tab20(np.arange(20)) if W.shape[1] > 10 else plt.cm.tab10(np.arange(10)); truth = conf.max() > 0
fig, axs = plt.subplots(1, 2 if truth else 1, figsize=(13 if truth else 6.5, 6)); axs = np.atleast_1d(axs)
for j in range(W.shape[1]): m = asg == j; axs[0].scatter(E[m, 0], E[m, 1], s=1, color=pal[j % len(pal)], label='state %d' % ids[j], rasterized=True)
axs[0].legend(markerscale=8, fontsize=8, frameon=False); axs[0].set_title('coloured by assigned state (cross-half memberships)', fontsize=10)
if truth:
    bc = np.bincount(conf[conf > 0]); cs = [k for k in np.argsort(-bc)[:8] if k > 0 and bc[k] > 0]
    m = conf == 0
    if m.any(): axs[1].scatter(E[m, 0], E[m, 1], s=1, color='#bbbbbb', label='junk', rasterized=True)
    m = (conf > 0) & ~np.isin(conf, cs); axs[1].scatter(E[m, 0], E[m, 1], s=1, color='#444444', label='other conformers', rasterized=True)
    for i, k in enumerate(cs): m = conf == k; axs[1].scatter(E[m, 0], E[m, 1], s=1, color=pal[i % 10], label='conformer %d' % k, rasterized=True)
    axs[1].legend(markerscale=8, fontsize=8, frameon=False); axs[1].set_title('coloured by the true conformer', fontsize=10)
for ax in axs: ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(a.title or 'per-particle UMAP of the cross-half state memberships', fontsize=11); plt.tight_layout(); plt.savefig(a.png, dpi=110, bbox_inches='tight'); print('png', a.png, len(W), 'particles')
