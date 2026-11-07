#!/usr/bin/env python3
"""Confusion-driven merge of a node's leaves from cross-half particle log-likelihoods (repro_r1/dscore_<run>; see decisive_merge.py), per run
and without reference to other runs: group score = max over its member leaves; each step merges the pair that is most often the top two of
an undecided particle (best minus second < THR nats), normalised by the pair's size; the partition is kept at every K (8 -> 2). Then the
runs' partitions compared at equal K (NMI, ARI, old-closed group Jaccard) next to the share of particles decided at that K.
usage: decisive_merge2.py RUN [RUN ...] [--thr 3]"""
import argparse, itertools, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='+'); ap.add_argument('--thr', type=float, default=3.); A = ap.parse_args()
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; H = '/mnt/beegfs/elmlund/afan/routes/toy55/simrnasep692'
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); oc = (old == 15) | (old == 16)
def nmi_ari(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((len(ua), len(ub))); np.add.at(C, (ia, ib), 1); n = C.sum(); P = C / n; pa, pb = P.sum(1), P.sum(0); nz = P > 0
    I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
    c2 = lambda x: x * (x - 1) / 2.; sij, sa, sb = c2(C).sum(), c2(C.sum(1)).sum(), c2(C.sum(0)).sum(); e = sa * sb / c2(n)
    return I / np.sqrt(max(h(pa) * h(pb), 1e-12)), (sij - e) / (0.5 * (sa + sb) - e)
def load(run):
    d = K_ + '/repro_r1/dscore_' + run; lab = np.loadtxt(d + '/labels.txt', dtype=int); S = np.full((len(lab), lab.max()), np.nan)
    for h in 'AB':
        mr = np.loadtxt(H + '_%s/meta_rows.txt' % h).astype(int) - 1; sc = np.load(d + '/scores_%s.npy' % h).astype(np.float64); m = lab[mr] > 0; S[mr[m]] = sc[m]
    return lab, S
def path(S, own, thr):
    groups = [[g] for g in range(S.shape[1])]; out = {}
    while True:
        G = np.stack([S[:, g].max(1) for g in groups], 1); o = np.argsort(-G, 1); t1, t2 = o[:, 0], o[:, 1]
        gap = G[np.arange(len(G)), t1] - G[np.arange(len(G)), t2]; gof = np.zeros(S.shape[1], int)
        for j, g in enumerate(groups):
            for leaf in g: gof[leaf] = j
        out[len(groups)] = (gof[own].copy(), (gap >= thr).mean())
        if len(groups) == 2: break
        und = gap < thr; cnt = np.zeros((len(groups), len(groups))); np.add.at(cnt, (np.minimum(t1[und], t2[und]), np.maximum(t1[und], t2[und])), 1)
        size = np.bincount(gof[own], minlength=len(groups)).astype(float); norm = cnt / np.maximum(np.add.outer(size, size), 1.)
        a, b = np.unravel_index(np.argmax(np.triu(norm, 1)), norm.shape); groups[a] = groups[a] + groups[b]; del groups[b]
    return out
D = {}; rows = None
for r in A.runs:
    lab, S = load(r); D[r] = (lab, S); ok = np.where(np.isfinite(S).all(1) & (lab > 0))[0]; rows = ok if rows is None else np.intersect1d(rows, ok)
P = {r: path(D[r][1][rows], D[r][0][rows] - 1, A.thr) for r in A.runs}; ocr = oc[rows]
print('%d runs, %d particles; confusion merge at %.1f nats' % (len(A.runs), len(rows), A.thr))
for k in range(8, 1, -1):
    line = 'K %d | decided %s |' % (k, ' '.join('%s %.0f%%' % (r, 100 * P[r][k][1]) for r in A.runs))
    for a, b in itertools.combinations(A.runs, 2):
        x, y = P[a][k][0], P[b][k][0]; nm, ar = nmi_ari(x, y)
        cx = np.unique(x)[np.argmax([ocr[x == g].mean() for g in np.unique(x)])]; cy = np.unique(y)[np.argmax([ocr[y == g].mean() for g in np.unique(y)])]
        X, Y = x == cx, y == cy; line += ' %s-%s NMI %.2f ARI %.2f closed J %.2f (%d/%d) |' % (a, b, nm, ar, (X & Y).sum() / (X | Y).sum(), X.sum(), Y.sum())
    print(line)
