#!/usr/bin/env python3
"""Site discovery, step 1: many small view-balanced groups of one node. The node's flex latent (all components) is centred per direction
bin (spherical k-means of the projection directions, ~PERBIN particles per bin; removes the view-dependent part of the latent) and
k-means-clustered into K = clamp(N / SIZE, 8, 99) groups, so a minority of a few % gets groups enriched in it (its site then shows in
the groups' cross-half variance even when the node's leading component is elsewhere). Writes groups.txt (per row, 0 = not a member).
Optional truth (--gt): the best concentration each conformer reaches in one group (share of the group, share of the conformer).
usage: site_groups.py --oris ptcl3D_dump.txt --rows rows.txt --flex flex_pca_coordinates.txt --out DIR [--size 8000 --perbin 1300 --gt GT]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--oris', required=True); ap.add_argument('--rows', required=True); ap.add_argument('--flex', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--size', type=float, default=8000.); ap.add_argument('--perbin', type=int, default=1300)
ap.add_argument('--k', type=int, default=0); ap.add_argument('--gt', default=''); ap.add_argument('--seed', type=int, default=0)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
rows = np.loadtxt(a.rows).astype(int); NR = len(rows); e1 = np.zeros(NR); e2 = np.zeros(NR); i = 0
for ln in open(a.oris):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); e1[i] = float(kv['e1']); e2[i] = float(kv['e2']); i += 1
assert i == NR, 'oris rows %d vs node rows %d' % (i, NR)
F = np.loadtxt(a.flex, comments='#', ndmin=2); fi = F[:, 0].astype(int); keep = rows[fi - 1] > 0; p = fi[keep]; Z = F[keep, 5:]; N = len(p)
E1, E2 = np.radians(e1[p - 1]), np.radians(e2[p - 1]); v = np.c_[np.sin(E2) * np.cos(E1), np.sin(E2) * np.sin(E1), np.cos(E2)]
nb = max(20, int(round(N / a.perbin))); C = v[rng.choice(N, nb, replace=False)].copy()
for it in range(40):
    lab = np.argmax(v @ C.T, 1)
    for k in range(nb):
        m = lab == k; C[k] = v[m].mean(0) / np.linalg.norm(v[m].mean(0)) if m.any() else v[rng.integers(N)]
lab = np.argmax(v @ C.T, 1); _, bi = np.unique(lab, return_inverse=True); nb = bi.max() + 1; cnt = np.bincount(bi)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=nb) for j in range(Z.shape[1])], 1) / cnt[:, None]; Zc = Z - bm[bi]; Zc -= Zc.mean(0)
Zc /= max(np.sqrt((Zc ** 2).sum(1).mean()), 1e-30)
K = a.k if a.k > 0 else int(np.clip(round(N / a.size), 8, 99))
best = None
for rep in range(3):   # k-means++ restarts, keep the lowest within-group sum of squares
    c0 = [Zc[rng.integers(N)]]
    for j in range(1, K):
        d2 = np.min([((Zc - c) ** 2).sum(1) for c in c0], 0); c0.append(Zc[rng.choice(N, p=d2 / d2.sum())])
    Cg = np.array(c0)
    for it in range(60):
        g = np.argmin(((Zc[:, None, :] - Cg[None]) ** 2).sum(2) if N * K < 4e7 else np.stack([((Zc - c) ** 2).sum(1) for c in Cg], 1), 1)
        for j in range(K): Cg[j] = Zc[g == j].mean(0) if (g == j).any() else Zc[rng.integers(N)]
    wss = sum(((Zc[g == j] - Cg[j]) ** 2).sum() for j in range(K))
    if best is None or wss < best[0]: best = (wss, g.copy())
g = best[1]; _, g = np.unique(g, return_inverse=True); K = g.max() + 1
out = np.zeros(NR, int); out[p - 1] = g + 1; np.savetxt(os.path.join(a.out, 'groups.txt'), out, fmt='%d')
print('node %d members, %d direction bins, %d-dim latent -> %d groups: sizes %s' % (N, nb, Z.shape[1], K, ' '.join(str(c) for c in np.bincount(g))))
if a.gt:
    G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(NR + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
    s = []
    for c in np.unique(conf):
        m = conf == c; sh = np.array([(m & (g == j)).sum() / max((g == j).sum(), 1) for j in range(K)]); j = int(sh.argmax())
        s.append((m.mean(), '%d: %.1f %% of the node -> best group %d %.0f %% of it (%.0f %% of the conformer)' % (c, 100 * m.mean(), j + 1, 100 * sh[j], 100 * (m & (g == j)).sum() / m.sum())))
    print('TRUTH per conformer (0 = junk), best concentration in one group:'); [print('   ' + t) for _, t in sorted(s, reverse=True)]
