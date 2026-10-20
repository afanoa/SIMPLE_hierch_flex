#!/usr/bin/env python3
"""Current-state figures. (A) full 692k RNase P run: UMAP of the responsibility vectors of an 80k subsample coloured by
density, z1, z2, soft-EM state, the sub-class closedness score (from closed_subclasses.py), tilt, junk mass, and the
sub-class mean positions. (B) toy10: per-conformer purity of the global latent vs the oracle ceiling vs base rate."""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S = '/private/tmp/claude-368467383/-Users-atifao-SIMPLE/25d47be4-7d1d-48d8-abdf-97b1c3c946ab/scratchpad'; rng = np.random.default_rng(0)
# ---- A
import umap
M = np.load(S + '/rnasep_klem/meta_all.npz'); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]; pos = {int(i): j for j, i in enumerate(ref)}
R = np.load(S + '/rnasep_klem/full_npz/klem_resp.npz'); r = R['r']; pinds = R['pinds'].astype(int); L = np.loadtxt(S + '/rnasep_klem/full_npz/klem_latent.txt', comments='#'); assert np.all(L[:, 0].astype(int) == pinds)
z1, z2, gm = L[:, 1], L[:, 2], L[:, 5]; n = len(pinds); sub = rng.choice(n, 80000, replace=False); rows = np.array([pos.get(int(p) - 1, -1) for p in pinds[sub]])
lab = np.where(rows >= 0, M['lab692'][np.maximum(rows, 0)], 0); st_all = np.loadtxt(S + '/rnasep_klem/states.txt').astype(int); st = np.where(rows >= 0, st_all[ref][np.maximum(rows, 0)], 0)
C = np.loadtxt('figs/fig_rnasep_subclass_closedness.txt'); cl = dict(zip(C[:, 0].astype(int), C[:, 4])); closed = np.array([cl.get(int(s), np.nan) for s in lab])
e2 = M['e2'][pinds[sub] - 1]
U = umap.UMAP(n_neighbors=30, min_dist=0.1, random_state=0).fit_transform(np.sqrt(r[sub].astype(float))); print('UMAP done')
fig, ax = plt.subplots(2, 4, figsize=(19, 9)); ax = ax.ravel(); kw = dict(s=1.2, linewidths=0, rasterized=True); g = gm[sub] > 0.5
ax[0].hexbin(U[:, 0], U[:, 1], gridsize=70, cmap='Greys', bins='log'); ax[0].set_title('full 692k run: responsibility UMAP (80k shown), density')
for k, (v, cm, t) in enumerate(((z1[sub], 'viridis', 'posterior-mean z1'), (z2[sub], 'viridis', 'posterior-mean z2'))):
    sc = ax[1 + k].scatter(U[g, 0], U[g, 1], c=v[g], cmap=cm, **kw); plt.colorbar(sc, ax=ax[1 + k]); ax[1 + k].set_title(t)
m = np.isfinite(closed); sc = ax[3].scatter(U[m, 0], U[m, 1], c=closed[m], cmap='coolwarm', vmin=-0.5, vmax=0.8, **kw); plt.colorbar(sc, ax=ax[3]); ax[3].set_title('closedness of the particle\'s 2D sub-class (vs deposited 18 - 28)')
for s_, c_ in zip((1, 2, 3), ('tab:green', 'tab:orange', 'tab:blue')):
    mm = st == s_; ax[4].scatter(U[mm, 0], U[mm, 1], color=c_, alpha=0.35, label='soft-EM state %d' % s_, **kw)
ax[4].legend(markerscale=8); ax[4].set_title('earlier soft-EM states')
sc = ax[5].scatter(U[:, 0], U[:, 1], c=e2, cmap='twilight', **kw); plt.colorbar(sc, ax=ax[5]); ax[5].set_title('tilt e2 (nuisance check)')
sc = ax[6].scatter(U[:, 0], U[:, 1], c=1 - gm[sub], cmap='magma', vmin=0, vmax=1, **kw); plt.colorbar(sc, ax=ax[6]); ax[6].set_title('mass in the free (junk) components')
# sub-class means on the grid, coloured by closedness, sized by members
A = np.loadtxt(S + '/rnasep_klem/full_npz/subclass_axis.txt'); su_, par_, n_, mt_, se_ = A.T
Lall = L; rows_all = np.array([pos.get(int(p) - 1, -1) for p in pinds]); lab_all = np.where(rows_all >= 0, M['lab692'][np.maximum(rows_all, 0)], 0); ok = (rows_all >= 0) & (gm > 0.5)
uu, ui = np.unique(lab_all[ok], return_inverse=True); cnt = np.bincount(ui); mz1 = np.bincount(ui, z1[ok]) / cnt; mz2 = np.bincount(ui, z2[ok]) / cnt; cc = np.array([cl.get(int(s), np.nan) for s in uu]); big = cnt >= 100
sc = ax[7].scatter(mz1[big], mz2[big], c=cc[big], cmap='coolwarm', vmin=-0.5, vmax=0.8, s=np.sqrt(cnt[big]) * 0.8, alpha=0.8, linewidths=0); plt.colorbar(sc, ax=ax[7]); ax[7].set_xlim(-0.5, 6.5); ax[7].set_ylim(-0.5, 6.5)
ax[7].set_title('sub-class MEAN positions on the 7x7 grid (size = members), closedness colour'); ax[7].set_xlabel('z1'); ax[7].set_ylabel('z2')
for x in ax[:7]: x.set_xticks([]); x.set_yticks([])
plt.tight_layout(); plt.savefig('figs/fig_state_rnasep_full_umap.png', dpi=105); print('figure -> figs/fig_state_rnasep_full_umap.png')
# ---- B: toy10 per-conformer purity: klem global (argmax one-to-one), oracle precision, base rate
from scipy.optimize import linear_sum_assignment
STATES = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; gt = {}
for l in open(S + '/toy10/gt_particles.txt'):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
def coverage(latent, nc):
    L = np.loadtxt(latent, comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; conf = np.array([gt.get(p, 0) for p in pind]); m = (gm > 0.5) & (conf > 0) & (am < nc)
    Cm = np.zeros((nc, len(STATES)))
    for j, c in enumerate(STATES): Cm[:, j] = np.bincount(am[m][conf[m] == c], minlength=nc)
    P = Cm / np.maximum(Cm.sum(1, keepdims=True), 1); tot = Cm.sum(1) > 0; r_, c_ = linear_sum_assignment(-P[tot]); pur = np.zeros(len(STATES)); pur[c_] = P[tot][r_, c_]
    best = P.max(0); return pur, best, Cm.sum(0) / Cm.sum()
pur, best, base = coverage(S + '/toy10/klem_latent.txt', 49)
O = np.load(S + '/toy10/oracle_lp10.txt.post.npz'); post, op = O['post'], O['pinds']; oc = np.array([gt.get(int(p), 0) for p in op]); om = oc > 0; pred = np.array(STATES)[post[om].argmax(1)]
oprec = np.array([np.mean(oc[om][pred == c] == c) if (pred == c).any() else 0. for c in STATES]); orec = np.array([np.mean(pred[oc[om] == c] == c) for c in STATES])
fig, ax = plt.subplots(1, 2, figsize=(15, 4.8)); x = np.arange(len(STATES)); wd = 0.27
ax[0].bar(x - wd, base, wd, color='lightgrey', label='base rate (chance)'); ax[0].bar(x, pur, wd, color='tab:blue', label='klem global latent, one cell per conformer (200k, 7x7)'); ax[0].bar(x + wd, oprec, wd, color='tab:red', label='oracle per-particle precision (truth maps, 10 A)')
ax[0].set_xticks(x); ax[0].set_xticklabels(STATES); ax[0].set_xlabel('toy10 conformer'); ax[0].set_ylabel('purity / precision'); ax[0].legend(fontsize=8); ax[0].set_title('every state has a cell; purity sits below the per-particle oracle ceiling'); ax[0].set_ylim(0, 1)
ax[1].bar(x - wd / 2, best / base, wd, color='tab:blue', label='klem: best-cell purity / base rate'); ax[1].bar(x + wd / 2, oprec / base, wd, color='tab:red', label='oracle precision / base rate')
ax[1].set_xticks(x); ax[1].set_xticklabels(STATES); ax[1].set_ylabel('enrichment over chance'); ax[1].legend(fontsize=8); ax[1].set_title('enrichment per conformer')
plt.tight_layout(); plt.savefig('figs/fig_state_toy10_coverage.png', dpi=110); print('figure -> figs/fig_state_toy10_coverage.png')
