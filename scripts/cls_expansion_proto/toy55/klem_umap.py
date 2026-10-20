#!/usr/bin/env python3
"""UMAP views of a klem run on REAL data (no truth): (1) UMAP of the per-particle responsibility vectors (49 grid cells +
free components, Hellinger metric) coloured by density, posterior-mean z1/z2, free-component (junk) mass, soft-EM
state, tilt, defocus and a second run's z1 (reproducibility); (2) UMAP of the flex_pca latent (z692) of the same
particles coloured by the klem z1 and by the soft-EM state (where the klem motion sits in the flex latent).
usage: klem_umap.py KLEMDIR META --labels states.txt [--other KLEMDIR2] [--out fig.png] [--out2 fig2.png]"""
import argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import umap
ap = argparse.ArgumentParser(); ap.add_argument('klemdir'); ap.add_argument('meta'); ap.add_argument('--labels', default=''); ap.add_argument('--other', default='')
ap.add_argument('--label_title', default='earlier soft-EM states (model output, not truth)'); ap.add_argument('--label_name', default='state'); ap.add_argument('--out', default='fig_klem_umap.png'); ap.add_argument('--out2', default='fig_flex_vs_klem_umap.png'); ap.add_argument('--nn', type=int, default=30); ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args(); M = np.load(a.meta); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]
R = np.load(a.klemdir + '/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); n = len(pinds); ngrid = len(R['i1'])
L = np.loadtxt(a.klemdir + '/klem_latent.txt', comments='#'); assert np.all(L[:, 0].astype(int) == pinds); z1, z2, gm = L[:, 1], L[:, 2], L[:, 5]
pind = pinds - 1; e2, dfx = M['e2'][pind], M['dfx'][pind]
pos = {int(i): j for j, i in enumerate(ref)}; rows = np.array([pos.get(int(i), -1) for i in pind]); ok = rows >= 0
st = None
if a.labels:
    lab = np.loadtxt(a.labels).astype(int); lab = lab[ref] if len(lab) == len(flags) else lab; st = np.where(ok, lab[np.maximum(rows, 0)], 0)
zo = None
if a.other:
    Lo = np.loadtxt(a.other + '/klem_latent.txt', comments='#'); d = dict(zip(Lo[:, 0].astype(int), zip(Lo[:, 1], Lo[:, 2], Lo[:, 5]))); t = np.array([d.get(int(i), (np.nan, np.nan, 0.)) for i in pinds]); zo, gmo = t[:, :2], t[:, 2]
print('%d particles; on-grid fraction %.3f; soft-EM labelled %.3f' % (n, np.mean(gm > 0.5), np.mean(st > 0) if st is not None else 0))
# 1. responsibility UMAP (Hellinger: euclidean on sqrt(r))
U = umap.UMAP(n_neighbors=a.nn, min_dist=0.1, metric='euclidean', random_state=a.seed).fit_transform(np.sqrt(r))
print('responsibility UMAP done'); np.savez_compressed(a.out.replace('.png', '_umap.npz'), U=U, pinds=pinds)
fig, ax = plt.subplots(2, 4, figsize=(19, 9)); ax = ax.ravel(); kw = dict(s=1.5, linewidths=0, rasterized=True)
ax[0].hexbin(U[:, 0], U[:, 1], gridsize=70, cmap='Greys', bins='log'); ax[0].set_title('responsibility UMAP, density (log)')
g = gm > 0.5
sc = ax[1].scatter(U[g, 0], U[g, 1], c=z1[g], cmap='viridis', **kw); plt.colorbar(sc, ax=ax[1]); ax[1].set_title('posterior-mean z1 (main latent axis)')
sc = ax[2].scatter(U[g, 0], U[g, 1], c=z2[g], cmap='viridis', **kw); plt.colorbar(sc, ax=ax[2]); ax[2].set_title('posterior-mean z2')
sc = ax[3].scatter(U[:, 0], U[:, 1], c=1 - gm, cmap='magma', vmin=0, vmax=1, **kw); plt.colorbar(sc, ax=ax[3]); ax[3].set_title('mass in the free (junk) components')
if st is not None:
    cmap_ = plt.get_cmap('tab10'); labs_ = sorted(set(st[st > 0]))
    for j_, s in enumerate(labs_):
        m = st == s; ax[4].scatter(U[m, 0], U[m, 1], color=cmap_(j_ % 10), alpha=0.35, label='%s %d (n=%d)' % (a.label_name, s, m.sum()), **kw); ax[4].text(np.median(U[m, 0]), np.median(U[m, 1]), str(s), fontsize=10, weight='bold')
    ax[4].legend(markerscale=8, fontsize=7); ax[4].set_title(a.label_title)
sc = ax[5].scatter(U[:, 0], U[:, 1], c=e2, cmap='twilight', **kw); plt.colorbar(sc, ax=ax[5]); ax[5].set_title('tilt e2 (nuisance check)')
sc = ax[6].scatter(U[:, 0], U[:, 1], c=dfx, cmap='cividis', **kw); plt.colorbar(sc, ax=ax[6]); ax[6].set_title('defocus (nuisance check)')
if zo is not None:
    # the grid has no preferred orientation across runs: compare the two 2-D latents by canonical correlation (orientation-free)
    # and the argmax cells by NMI; colour by the other run's first canonical variate
    m = g & np.isfinite(zo[:, 0]) & (gmo > 0.5)
    X = np.c_[z1, z2][m]; Y = zo[m]; X = X - X.mean(0); Y = Y - Y.mean(0)
    Sxx, Syy, Sxy = X.T @ X, Y.T @ Y, X.T @ Y; ix, iy = np.linalg.inv(Sxx), np.linalg.inv(Syy)
    w, v = np.linalg.eig(ix @ Sxy @ iy @ Sxy.T); o = np.argsort(-w.real); cc = np.sqrt(np.clip(w.real[o], 0, 1)); ax_ = v[:, o[0]].real
    b = iy @ Sxy.T @ ax_; yv = Y @ b; yv = (yv - yv.min()) / (yv.max() - yv.min())
    def nmi(x, y):
        x = np.unique(x, return_inverse=True)[1]; y = np.unique(y, return_inverse=True)[1]; C = np.zeros((x.max() + 1, y.max() + 1)); np.add.at(C, (x, y), 1); P = C / C.sum(); px, py = P.sum(1), P.sum(0)
        with np.errstate(divide='ignore', invalid='ignore'): I = np.nansum(P * np.log(P / (px[:, None] * py[None, :])))
        hx = -np.sum(px[px > 0] * np.log(px[px > 0])); hy = -np.sum(py[py > 0] * np.log(py[py > 0])); return I / max(np.sqrt(hx * hy), 1e-12)
    cells = np.clip(np.rint(z1), 0, 6) * 7 + np.clip(np.rint(z2), 0, 6); cello = np.clip(np.rint(zo[:, 0]), 0, 6) * 7 + np.clip(np.rint(zo[:, 1]), 0, 6)
    print('two runs: canonical correlations %.2f %.2f; NMI(cells) %.3f; junk agreement %.3f' % (cc[0], cc[1], nmi(cells[m], cello[m]), np.mean((gm > 0.5) == (gmo > 0.5))))
    sc = ax[7].scatter(U[m, 0], U[m, 1], c=yv, cmap='viridis', **kw); plt.colorbar(sc, ax=ax[7])
    ax[7].set_title('RANDOM-start run, 1st canonical variate (canon. corr %.2f, %.2f)' % (cc[0], cc[1]))
for x in ax: x.set_xticks([]); x.set_yticks([])
plt.tight_layout(); plt.savefig(a.out, dpi=110); print('figure ->', a.out)
# 2. flex latent UMAP of the same particles
if 'z692' in M.files and ok.sum() > 1000:
    Z = M['z692'][rows[ok]]; Z = (Z - Z.mean(0)) / Z.std(0)
    V = umap.UMAP(n_neighbors=a.nn, min_dist=0.1, random_state=a.seed).fit_transform(Z); print('flex latent UMAP done')
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
    ax[0].hexbin(V[:, 0], V[:, 1], gridsize=70, cmap='Greys', bins='log'); ax[0].set_title('flex_pca latent (16-D) UMAP, density')
    gg = g[ok]; sc = ax[1].scatter(V[gg, 0], V[gg, 1], c=z1[ok][gg], cmap='viridis', **kw); plt.colorbar(sc, ax=ax[1]); ax[1].set_title('coloured by klem z1')
    sc = ax[2].scatter(V[:, 0], V[:, 1], c=(1 - gm)[ok], cmap='magma', vmin=0, vmax=1, **kw); plt.colorbar(sc, ax=ax[2]); ax[2].set_title('coloured by klem junk mass')
    if st is not None:
        so = st[ok]
        for s, c in zip(sorted(set(so[so > 0])), ('tab:green', 'tab:orange', 'tab:blue', 'tab:red')):
            m = so == s; ax[3].scatter(V[m, 0], V[m, 1], color=c, alpha=0.35, label='state %d' % s, **kw)
        ax[3].legend(markerscale=8); ax[3].set_title('soft-EM states')
    for x in ax: x.set_xticks([]); x.set_yticks([])
    plt.tight_layout(); plt.savefig(a.out2, dpi=110); print('figure ->', a.out2)
