#!/usr/bin/env python3
"""Particle-level width rule for a tree node, tested across independent runs: per run, cross-half log-likelihoods of every node particle
against that run's K leaf maps (score_embed.py: repro_r1/dscore_<run>/scores_{A,B}.npy, labels.txt = the run's leaf per master row);
leaves merged greedily, least separated pair first, while the pair's median per-particle margin (own group vs the other, group score = max
over its member leaves) is below THR nats; then the merged partitions compared across runs (NMI, ARI, old-closed group Jaccard), with the
share of particles decided by >= THR nats per run.  usage: decisive_merge.py RUN [RUN ...] [--thr 2,3,4]"""
import sys, itertools, numpy as np
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; H = '/mnt/beegfs/elmlund/afan/routes/toy55/simrnasep692'
import argparse
ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='+'); ap.add_argument('--thr', default='2,3,4'); A = ap.parse_args()
runs = A.runs; THR = [float(x) for x in A.thr.split(',')]
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); oc = (old == 15) | (old == 16)
def nmi_ari(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((len(ua), len(ub))); np.add.at(C, (ia, ib), 1); n = C.sum(); P = C / n; pa, pb = P.sum(1), P.sum(0); nz = P > 0
    I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
    c2 = lambda x: x * (x - 1) / 2.; sij, sa, sb = c2(C).sum(), c2(C.sum(1)).sum(), c2(C.sum(0)).sum(); e = sa * sb / c2(n)
    return I / np.sqrt(max(h(pa) * h(pb), 1e-12)), (sij - e) / (0.5 * (sa + sb) - e)
def load(run):
    d = K_ + '/repro_r1/dscore_' + run; lab = np.loadtxt(d + '/labels.txt', dtype=int); k = lab.max(); S = np.full((len(lab), k), np.nan)
    for h in 'AB':
        mr = np.loadtxt(H + '_%s/meta_rows.txt' % h).astype(int) - 1; sc = np.load(d + '/scores_%s.npy' % h).astype(np.float64); m = lab[mr] > 0; S[mr[m]] = sc[m]
    rows = np.where(np.isfinite(S).all(1) & (lab > 0))[0]; return lab, S, rows
def merge(S, own, thr):
    groups = [[g] for g in range(S.shape[1])]; gof = np.arange(S.shape[1]); hist = []
    while len(groups) > 1:
        G = np.stack([S[:, g].max(1) for g in groups], 1); cur = gof[own]; best = None
        for a, b in itertools.combinations(range(len(groups)), 2):
            ma, mb = cur == a, cur == b
            sep = np.median(np.r_[G[ma, a] - G[ma, b], G[mb, b] - G[mb, a]])
            if best is None or sep < best[0]: best = (sep, a, b)
        if best[0] >= thr: break
        _, a, b = best; hist.append('%.1f' % best[0]); groups[a] = groups[a] + groups[b]; del groups[b]
        for j, g in enumerate(groups):
            for leaf in g: gof[leaf] = j
    G = np.stack([S[:, g].max(1) for g in groups], 1); srt = np.sort(G, 1); dec = (srt[:, -1] - srt[:, -2] >= thr).mean() if G.shape[1] > 1 else 1.
    return gof[own], groups, dec, hist
D = {r: load(r) for r in runs}; common = None
for r in runs: common = D[r][2] if common is None else np.intersect1d(common, D[r][2])
print('%d runs, %d particles scored in all' % (len(runs), len(common)))
for thr in THR:
    print('== threshold %.1f nats (merge while the least separated pair has a median margin below it)' % thr); ML = {}
    for r in runs:
        lab, S, rows = D[r]; own = lab[common] - 1; m, groups, dec, hist = merge(S[common], own, thr); ML[r] = m
        sizes = np.bincount(m); en = [oc[common][m == g].mean() / oc[common].mean() for g in range(len(groups))]
        print('  %-8s %d -> %d groups, sizes %s, closed enrichment %s, decided >= %.0f nats %.0f %%; merge margins %s' % (r, S.shape[1], len(groups),
              ' '.join(str(x) for x in sizes), ' '.join('%.1f' % x for x in en), thr, 100 * dec, ' '.join(hist)))
    for a, b in itertools.combinations(runs, 2):
        x, y = ML[a], ML[b]; nm, ar = nmi_ari(x, y)
        ca = np.argmax([oc[common][x == g].mean() for g in np.unique(x)]); cb = np.argmax([oc[common][y == g].mean() for g in np.unique(y)])
        A_, B_ = set(np.where(x == np.unique(x)[ca])[0]), set(np.where(y == np.unique(y)[cb])[0])
        print('  %s vs %s: NMI %.2f ARI %.2f, closed groups %d / %d Jaccard %.2f' % (a, b, nm, ar, len(A_), len(B_), len(A_ & B_) / len(A_ | B_)))
