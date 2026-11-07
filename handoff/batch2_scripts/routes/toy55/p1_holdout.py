#!/usr/bin/env python3
"""Held-out comparison of a node's leaf models (problem 1) from the cross-half scoring passes. Per particle and model the best distance
over that model's NB maps (best-candidate euclid dist of the refine3D prob tables; even particles from pass E = odd-half maps, odd ones
from pass O = even-half maps, so no particle meets a map that contains it). Lower = the model explains the particle better. Reports the
gain of the hard-state and coordinate-bin models over the random split (no heterogeneity), their paired difference with its standard
error, per half; with truth: NMI(conformer), swing eta2 of each partition and of the held-out best-map assignment. Also the score-PCA
coordinate (lrec_pc2: centred scores, PC1-2) from the cross-half hard block against the self-inclusive pass S and the old it0 PCs.
usage: p1_holdout.py OUT --passE DIR --passO DIR [--passS DIR] --nb 8 [--gt GT --it0 coord_dir/it0/pc]"""
import os, glob, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--passE', required=True); ap.add_argument('--passO', required=True); ap.add_argument('--passS', default='')
ap.add_argument('--nb', type=int, default=8); ap.add_argument('--gt', default=''); ap.add_argument('--it0', default=''); a = ap.parse_args()
CAND = np.dtype([('iref', '<i4'), ('inpl', '<i4'), ('dist', '<f4'), ('x', '<f4'), ('y', '<f4'), ('has_sh', '<i4')])
def tables(d, nst, nrows):   # per-row best dist to every state (inf = not scored), as lrec_pc2.py reads them
    D = np.full((nrows, nst), np.inf)
    for f in sorted(glob.glob('%s/*_refine3D/dist_part_neigh_*.dat' % d)):
        sz = os.path.getsize(f); raw = open(f, 'rb').read(); nrefs, nptcls, nnz, nchunks = np.frombuffer(raw, '<i8', 4); nspace = nrefs // nst; assert nspace * nst == nrefs
        pinds = np.frombuffer(raw, '<i4', nptcls, 32); pos = sz - (28 * nnz + 4 * nchunks); Dp = np.full((nptcls, nst), np.inf)
        for c in range(nchunks):
            m = int(np.frombuffer(raw, '<i4', 1, pos)[0]); pos += 4; idx = np.frombuffer(raw, '<i4', m, pos) - 1; pos += 4 * m
            cand = np.frombuffer(raw, CAND, m, pos); pos += 24 * m; st = (cand['iref'] - 1) // nspace; ok = (cand['inpl'] > 0) & (st >= 0) & (st < nst)
            np.minimum.at(Dp, (idx[ok], st[ok]), cand['dist'][ok])
        assert pos == sz, (f, pos, sz); D[pinds - 1] = Dp
    return D
P = os.path.dirname(os.path.abspath(a.passE)); ld = lambda nm: np.loadtxt('%s/%s.txt' % (P, nm), dtype=int)
hard, bins, rnd, eo, com = ld('hard'), ld('bins'), ld('rand'), ld('eo'), ld('common') > 0; nr = len(hard); nb = a.nb
DE, DO = tables(a.passE, 3 * nb, nr), tables(a.passO, 3 * nb, nr); D = np.where((eo == 0)[:, None], DE, DO); ok = com & np.isfinite(D).all(1)
h, b, r = D[:, :nb].min(1), D[:, nb:2 * nb].min(1), D[:, 2 * nb:].min(1)
out = ['%d common rows, %d fully scored (even %d, odd %d); distances: lower = better fit' % (com.sum(), ok.sum(), (ok & (eo == 0)).sum(), (ok & (eo == 1)).sum())]
def line(tag, m):
    gh, gb, dd = (r - h)[m], (r - b)[m], (b - h)[m]; se = dd.std() / np.sqrt(m.sum())
    return '%-5s n %6d | gain over random: hard %.4f  bins %.4f | bins - hard %+.4f +- %.4f (z %+.1f) | hard better for %.1f %% of particles' % (tag, m.sum(), gh.mean(), gb.mean(), dd.mean(), se, dd.mean() / max(se, 1e-300), 100. * (h[m] < b[m]).mean())
out += [line('all', ok), line('even', ok & (eo == 0)), line('odd', ok & (eo == 1))]
G = np.loadtxt(a.gt, comments='#') if a.gt else None
if G is not None:
    conf = np.full(nr, -1); gi = G[:, 0].astype(int) - 1; v = (gi >= 0) & (gi < nr); conf[gi[v]] = G[v, 1].astype(int); sw = ok & (conf >= 200) & (conf <= 230); x = (conf - 200.)
    def nmi(u, w):
        _, ui = np.unique(u, return_inverse=True); _, wi = np.unique(w, return_inverse=True); J = np.zeros((ui.max() + 1, wi.max() + 1)); np.add.at(J, (ui, wi), 1); J /= J.sum()
        pu, pw = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(pu, pw)[nz])).sum(); return I / max(np.sqrt((pu * np.log(pu)).sum() * (pw * np.log(pw)).sum()), 1e-12)
    def eta(l):
        xs, ls = x[sw], l[sw]
        if sw.sum() < 50: return np.nan
        mu = xs.mean(); return sum((ls == s).sum() * (xs[ls == s].mean() - mu) ** 2 for s in np.unique(ls)) / max(((xs - mu) ** 2).sum(), 1e-12)
    ah, ab = D[:, :nb].argmin(1) + 1, D[:, nb:2 * nb].argmin(1) + 1
    out.append('TRUTH (%d swing particles = %.0f %%): partition NMI(conformer) / swing eta2: hard %.3f / %.3f, bins %.3f / %.3f, rand %.3f / %.3f; held-out best-map assignment: hard %.3f / %.3f, bins %.3f / %.3f'
               % (sw.sum(), 100. * sw.sum() / ok.sum(), nmi(conf[ok], hard[ok]), eta(hard), nmi(conf[ok], bins[ok]), eta(bins), nmi(conf[ok], rnd[ok]), eta(rnd), nmi(conf[ok], ah[ok]), eta(ah), nmi(conf[ok], ab[ok]), eta(ab)))
    def pcs(Dh, m):   # lrec_pc2's coordinate: per-particle centred, then column-centred scores -> PC1, PC2
        X = Dh[m] - Dh[m].mean(1, keepdims=True); X -= X.mean(0); U, S, Vt = np.linalg.svd(X, full_matrices=False); Z = np.full((nr, 2), np.nan); Z[m] = U[:, :2] * S[:2]; return Z, S[:5] ** 2 / (S ** 2).sum()
    def cor(Z):
        m = sw & np.isfinite(Z).all(1); wv = np.linalg.lstsq(np.c_[Z[m], np.ones(m.sum())], x[m], rcond=None)[0][:2]
        return '|corr| PC1 %.3f PC2 %.3f, best direction %.3f' % (abs(np.corrcoef(Z[m, 0], x[m])[0, 1]), abs(np.corrcoef(Z[m, 1], x[m])[0, 1]), np.corrcoef(Z[m] @ wv, x[m])[0, 1])
    Zx, fx = pcs(D[:, :nb], ok); out.append('score coordinate vs swing (swing particles): cross-half hard maps: %s (variance shares %s)' % (cor(Zx), ' '.join('%.3f' % q for q in fx[:3])))
    if a.passS:
        DS = tables(a.passS, nb, nr); mS = com & np.isfinite(DS).all(1); ZS, fS = pcs(DS, mS)
        out.append('                                          self-inclusive (same maps, averages): %s (variance shares %s)' % (cor(ZS), ' '.join('%.3f' % q for q in fS[:3])))
    if a.it0 and os.path.isfile(a.it0 + '/pc1.txt'):
        Z0 = np.c_[np.loadtxt(a.it0 + '/pc1.txt'), np.loadtxt(a.it0 + '/pc2.txt')]; out.append('                                          old it0 (self-inclusive, nspace 20000): %s' % cor(Z0))
open(a.out + '/holdout.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
