#!/usr/bin/env python3
"""p1_holdout.py's tables re-read for the 1-D capture ratio (problem 1): rho = gain(bins) / gain(hard) over the random split = the share of
the heterogeneity the hard states explain on held-out particles that NB equal-width bins of the node's 1-D coordinate also explain
(1 = the node is one axis; small = several independent directions or identities), with a particle bootstrap; and per truth class
(swing 200-230 / other conformers / junk) the two gains (where each model wins).
usage: p1_holdout2.py NODEDIR [NODEDIR ...] --nb 8 [--gt GT --nboot 200]"""
import os, glob, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('dirs', nargs='+'); ap.add_argument('--nb', type=int, default=8); ap.add_argument('--gt', default=''); ap.add_argument('--nboot', type=int, default=200)
a = ap.parse_args(); CAND = np.dtype([('iref', '<i4'), ('inpl', '<i4'), ('dist', '<f4'), ('x', '<f4'), ('y', '<f4'), ('has_sh', '<i4')])
def tables(d, nst, nrows):
    D = np.full((nrows, nst), np.inf)
    for f in sorted(glob.glob('%s/*_refine3D/dist_part_neigh_*.dat' % d)):
        sz = os.path.getsize(f); raw = open(f, 'rb').read(); nrefs, nptcls, nnz, nchunks = np.frombuffer(raw, '<i8', 4); nspace = nrefs // nst
        pinds = np.frombuffer(raw, '<i4', nptcls, 32); pos = sz - (28 * nnz + 4 * nchunks); Dp = np.full((nptcls, nst), np.inf)
        for c in range(nchunks):
            m = int(np.frombuffer(raw, '<i4', 1, pos)[0]); pos += 4; idx = np.frombuffer(raw, '<i4', m, pos) - 1; pos += 4 * m
            cand = np.frombuffer(raw, CAND, m, pos); pos += 24 * m; st = (cand['iref'] - 1) // nspace; ok = (cand['inpl'] > 0) & (st >= 0) & (st < nst); np.minimum.at(Dp, (idx[ok], st[ok]), cand['dist'][ok])
        D[pinds - 1] = Dp
    return D
G = np.loadtxt(a.gt, comments='#') if a.gt else None; rng = np.random.default_rng(1); nb = a.nb
for d in a.dirs:
    ld = lambda nm: np.loadtxt('%s/%s.txt' % (d, nm), dtype=int); eo, com = ld('eo'), ld('common') > 0; nr = len(eo)
    D = np.where((eo == 0)[:, None], tables(d + '/pass_E', 3 * nb, nr), tables(d + '/pass_O', 3 * nb, nr)); ok = com & np.isfinite(D).all(1)
    gh, gb = (D[:, 2 * nb:].min(1) - D[:, :nb].min(1))[ok], (D[:, 2 * nb:].min(1) - D[:, nb:2 * nb].min(1))[ok]
    bo = [gb[i].mean() / gh[i].mean() for i in (rng.integers(0, len(gh), len(gh)) for _ in range(a.nboot))]
    s = '%s: %d particles | gain hard %.6f bins %.6f | rho = bins / hard %.3f +- %.3f' % (os.path.basename(d), ok.sum(), gh.mean(), gb.mean(), gb.mean() / gh.mean(), np.std(bo))
    if G is not None:
        conf = np.full(nr, -1); gi = G[:, 0].astype(int) - 1; v = (gi >= 0) & (gi < nr); conf[gi[v]] = G[v, 1].astype(int); c = conf[ok]
        for nm, m in (('swing', (c >= 200) & (c <= 230)), ('other', (c > 0) & ((c < 200) | (c > 230))), ('junk', c == 0)):
            if m.sum() >= 100: s += ' | %s %d: hard %.6f bins %.6f (rho %.2f)' % (nm, m.sum(), gh[m].mean(), gb[m].mean(), gb[m].mean() / max(gh[m].mean(), 1e-12))
    print(s, flush=True)
