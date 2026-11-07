#!/usr/bin/env python3
"""Sanity readouts of a klem latent on REAL data (no truth): where the earlier soft-EM states (S1/S2/S4), the
cls_expansion sub-classes, the VIEW (tilt/rotation bins) and the DEFOCUS sit on the latent grid. A latent that
follows view or defocus is a nuisance latent; one that follows the states / sub-classes but not view is structural.
usage: klem_real_readout.py KLEMDIR META [--labels states.txt] [--n1 7 --n2 7] [--out fig.png]"""
import sys, os, argparse, numpy as np
from scipy.stats import spearmanr
ap = argparse.ArgumentParser(); ap.add_argument('klemdir'); ap.add_argument('meta'); ap.add_argument('--labels', default=''); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); ap.add_argument('--out', default='')
a = ap.parse_args(); M = np.load(a.meta); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]
L = np.loadtxt(a.klemdir + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int) - 1; z1, z2, am, maxr, gm = L[:, 1], L[:, 2], L[:, 3].astype(int), L[:, 4], L[:, 5]
ngrid = a.n1 * a.n2; ongrid = gm > 0.5; print('%d particles, %.3f on the grid (rest in free components), mean max r %.2f' % (len(pind), ongrid.mean(), maxr.mean()))
def nmi(x, y):
    x = np.unique(x, return_inverse=True)[1]; y = np.unique(y, return_inverse=True)[1]; C = np.zeros((x.max() + 1, y.max() + 1)); np.add.at(C, (x, y), 1); P = C / C.sum(); px, py = P.sum(1), P.sum(0)
    with np.errstate(divide='ignore', invalid='ignore'): I = np.nansum(P * np.log(P / (px[:, None] * py[None, :])))
    hx = -np.sum(px[px > 0] * np.log(px[px > 0])); hy = -np.sum(py[py > 0] * np.log(py[py > 0])); return I / max(np.sqrt(hx * hy), 1e-12)
def cell_of(z1v, z2v): return np.clip(np.rint(z1v), 0, a.n1 - 1).astype(int) * a.n2 + np.clip(np.rint(z2v), 0, a.n2 - 1).astype(int)
cells = cell_of(z1, z2); g = ongrid
# nuisance: view (e2 tilt bins and e1 bins), defocus
e1, e2, dfx = M['e1'][pind], M['e2'][pind], M['dfx'][pind]
vb = np.digitize(e2, np.percentile(e2, np.linspace(0, 100, 9)[1:-1])); rb = np.digitize(e1, np.linspace(0, 360, 9)[1:-1]); db = np.digitize(dfx, np.percentile(dfx, np.linspace(0, 100, 9)[1:-1]))
print('NMI(cell, tilt-bin) %.3f | NMI(cell, rotation-bin) %.3f | NMI(cell, defocus-bin) %.3f | Spearman(z1, tilt) %+.2f (z2 %+.2f) | Spearman(z1, defocus) %+.2f (z2 %+.2f)' %
      (nmi(cells[g], vb[g]), nmi(cells[g], rb[g]), nmi(cells[g], db[g]), spearmanr(z1[g], e2[g])[0], spearmanr(z2[g], e2[g])[0], spearmanr(z1[g], dfx[g])[0], spearmanr(z2[g], dfx[g])[0]))
# sub-classes (lab692 rows = refined subset order)
if 'lab692' in M.files:
    pos = {int(i): j for j, i in enumerate(ref)}; rows = np.array([pos.get(int(i), -1) for i in pind]); ok = rows >= 0; sub = M['lab692'][rows[ok]]; par = M['parent692'][rows[ok]]
    rng = np.random.default_rng(0); subp = sub.copy()   # null: sub-class labels shuffled WITHIN each parent (view structure kept, within-view structure destroyed)
    for p in np.unique(par): m = par == p; subp[m] = rng.permutation(sub[m])
    print('NMI(cell, sub-class) %.3f (within-parent shuffle null %.3f) | NMI(cell, parent class) %.3f (shuffle null %.3f)  (sub-classes share view with their parent: the difference is the within-class part; NMI has a finite-sample bias, compare with the nulls)' %
          (nmi(cells[ok & g], sub[g[ok]]), nmi(cells[ok & g], subp[g[ok]]), nmi(cells[ok & g], par[g[ok]]), nmi(cells[ok & g], rng.permutation(par[g[ok]]))))
    # within-parent: does the latent separate sub-classes of the SAME parent (view held fixed)? with the within-parent shuffle null
    vals = []; nulls = []
    for p in np.unique(par):
        m = (par == p) & g[ok]; s = sub[m]
        if len(np.unique(s)) > 2 and m.sum() > 200: c_ = cells[ok & g][par[g[ok]] == p]; vals.append(nmi(c_, s)); nulls.append(nmi(c_, rng.permutation(s)))
    print('within-parent NMI(cell, sub-class): median %.3f over %d parents (view fixed); within-parent shuffle null median %.3f' % (np.median(vals), len(vals), np.median(nulls)))
if a.labels:
    lab = np.loadtxt(a.labels).astype(int); lab = lab[ref] if len(lab) == len(flags) else lab; st = lab[rows[ok]] if 'lab692' in M.files else None
    if st is not None:
        m = (st > 0) & g[ok]; print('soft-EM states: NMI(cell, state) %.3f (shuffle null %.3f); mean (z1, z2) per state: %s' % (nmi(cells[ok & g][st[g[ok]] > 0], st[m]), nmi(cells[ok & g][st[g[ok]] > 0], np.random.default_rng(1).permutation(st[m])), ' '.join('%d:(%.2f,%.2f) n=%d' % (s, z1[ok & g][st[g[ok]] == s].mean(), z2[ok & g][st[g[ok]] == s].mean(), np.sum(st[g[ok]] == s)) for s in np.unique(st[m]))))
if a.out:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    H, _, _ = np.histogram2d(z1[g], z2[g], bins=[np.linspace(-0.5, a.n1 - 0.5, 4 * a.n1 + 1), np.linspace(-0.5, a.n2 - 0.5, 4 * a.n2 + 1)]); ax[0].imshow(H.T, origin='lower', extent=(-0.5, a.n1 - 0.5, -0.5, a.n2 - 0.5), cmap='Greys', aspect='auto'); ax[0].set_title('posterior-mean positions')
    sc = ax[1].scatter(z1[g][:20000], z2[g][:20000], c=e2[g][:20000], s=2, cmap='viridis'); plt.colorbar(sc, ax=ax[1]); ax[1].set_title('coloured by tilt e2')
    if a.labels and st is not None:
        for s in np.unique(st[m]): mm = (st == s) & g[ok]; ax[2].scatter(z1[ok][mm][:5000], z2[ok][mm][:5000], s=2, alpha=0.3, label='state %d' % s)
        ax[2].legend(); ax[2].set_title('soft-EM states')
    for x in ax: x.set_xlim(-0.5, a.n1 - 0.5); x.set_ylim(-0.5, a.n2 - 0.5)
    plt.tight_layout(); plt.savefig(a.out, dpi=110); print('figure ->', a.out)
