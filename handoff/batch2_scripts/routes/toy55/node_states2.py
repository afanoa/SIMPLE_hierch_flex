#!/usr/bin/env python3
"""States of a klem board without any truth: surviving cells (cross-half reproducibility >= --minrep) clustered by the
cross-half similarity of their deviations (average linkage at --sim); every particle is assigned to the cluster that
holds most of its responsibility mass and, within each cluster, the top N (N = min(--nmax, --frac x cluster mass)) by
that mass are selected; writes OUT/states_clusters.txt (one label per project row, 0 = unselected), OUT/clusters.txt and
OUT/cluster_targets.npz (each cluster's mean deviation on the band, the level-2 ranking target).
--support: a cluster is a state only if the data support it: its mean deviation reproduces between the halves at >= --sim
and it is separable (cross-half similarity to every other cluster < --sim); no particle count enters. A cluster that fails
is merged into its most similar neighbour if the merged cluster passes, otherwise it writes no state.
usage: node_states2.py BOARD OUT --nrows N [--sim 0.7 --nmax 30000 --frac 0.8 --minmass 2000 --maxk 8 --support]"""
import sys, os, argparse, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('out'); ap.add_argument('--nrows', type=int, required=True); ap.add_argument('--sim', type=float, default=0.7)
ap.add_argument('--nmax', type=int, default=30000); ap.add_argument('--frac', type=float, default=0.8); ap.add_argument('--minmass', type=float, default=2000.); ap.add_argument('--maxk', type=int, default=8); ap.add_argument('--minrep', type=float, default=0.5); ap.add_argument('--support', action='store_true'); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); H = np.load(a.board + '/klem_halves.npz'); D0, D1 = H['D0'], H['D1']; R = np.load(a.board + '/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(D0); a.nrows = a.nrows or int(pinds.max())
st = np.loadtxt(a.board + '/klem_states.txt', comments='#'); crep = st[:ngrid, 2]; keep = np.where(crep >= a.minrep)[0]; n = len(keep)
S = np.zeros((n, n))
for ii, i in enumerate(keep):
    for jj, j in enumerate(keep):
        c01 = np.real(np.sum(np.conj(D0[i]) * D1[j])) / max(np.sqrt(np.sum(np.abs(D0[i]) ** 2) * np.sum(np.abs(D1[j]) ** 2)), 1e-30); c10 = np.real(np.sum(np.conj(D1[i]) * D0[j])) / max(np.sqrt(np.sum(np.abs(D1[i]) ** 2) * np.sum(np.abs(D0[j]) ** 2)), 1e-30)
        S[ii, jj] = 0.5 * (c01 + c10)
Dm = np.clip(1 - S, 0, 2); np.fill_diagonal(Dm, 0); Z = linkage(squareform(Dm, checks=False), method='average'); lab = fcluster(Z, t=1 - a.sim, criterion='distance'); K = lab.max()
mass_c = np.array([r[:, keep[lab == c]].sum() for c in range(1, K + 1)]); order = np.argsort(-mass_c)
def xh(cells_a, cells_b):   # cross-half similarity of two cell groups' mean deviations (symmetrised)
    a0, a1, b0, b1 = D0[cells_a].mean(0), D1[cells_a].mean(0), D0[cells_b].mean(0), D1[cells_b].mean(0)
    cc = lambda x, y: np.real(np.sum(np.conj(x) * y)) / max(np.sqrt(np.sum(np.abs(x) ** 2) * np.sum(np.abs(y) ** 2)), 1e-30)
    return 0.5 * (cc(a0, b1) + cc(a1, b0))
def supported(c):
    cells = keep[lab == c]; rep = xh(cells, cells); sep = max([xh(cells, keep[lab == d]) for d in range(1, K + 1) if d != c] or [-1.])
    return rep >= a.sim and sep < a.sim, rep, sep
dropped = []
while a.support and K > 1:
    bad = [(c, ) + supported(c) for c in range(1, K + 1)]; bad = [b for b in bad if not b[1]]
    if not bad: break
    c, _, rep, sep = min(bad, key=lambda b: mass_c[b[0] - 1]); others = [d for d in range(1, K + 1) if d != c]
    d = max(others, key=lambda d: xh(keep[lab == c], keep[lab == d])); lab2 = lab.copy(); lab2[lab2 == c] = d
    m_rep = xh(keep[lab2 == d], keep[lab2 == d]); m_sep = max([xh(keep[lab2 == d], keep[lab2 == e]) for e in range(1, K + 1) if e not in (c, d)] or [-1.])
    if m_rep >= a.sim and m_sep < a.sim:
        print('  cluster of %d cells (mass %.0f) unsupported: reproducibility %.2f separability %.2f -> merged into its neighbour (%d cells, mass %.0f); merged: %.2f / %.2f' % ((lab == c).sum(), mass_c[c - 1], rep, sep, (lab == d).sum(), mass_c[d - 1], m_rep, m_sep))
        lab[lab == c] = d
    else:
        print('  cluster of %d cells (mass %.0f) unsupported: reproducibility %.2f separability %.2f -> no state (merge would give %.2f / %.2f)' % ((lab == c).sum(), mass_c[c - 1], rep, sep, m_rep, m_sep)); dropped.append(keep[lab == c]); keep_mask = lab != c; lab = lab[keep_mask]; keep = keep[keep_mask]
    lab[lab > c] -= 1; K -= 1
    mass_c = np.array([r[:, keep[lab == e]].sum() for e in range(1, K + 1)])
order = np.argsort(-mass_c)
print('board %s: %d of %d cells survive, %d clusters at similarity %.1f; cluster masses (particles): %s' % (a.board.rstrip('/').split('/')[-1], n, ngrid, K, a.sim, ' '.join('%d' % mass_c[c] for c in order[:12])))
M = np.stack([r[:, keep[lab == c]].sum(1) for c in range(1, K + 1)], 1); am = M.argmax(1); states = np.zeros(a.nrows, int); out = []; T = []
for rank, c in enumerate(order[:a.maxk]):
    if mass_c[c] < a.minmass: continue
    mine = np.where(am == c)[0]; N = int(min(a.nmax, a.frac * mass_c[c], len(mine))); top = mine[np.argsort(-M[mine, c])[:N]]; states[pinds[top] - 1] = rank + 1
    cells = keep[lab == c + 1]; T.append(np.mean(0.5 * (D0[cells] + D1[cells]), 0)); out.append((rank + 1, len(cells), mass_c[c], N, ' '.join('%d' % (k + 1) for k in cells)))
    print('  state %d: %d cells, mass %.0f, selected %d particles (mass in the cluster of the selected: median %.2f)' % (rank + 1, len(cells), mass_c[c], N, np.median(M[top, c])))
np.savetxt(a.out + '/states_clusters.txt', states, fmt='%d')
with open(a.out + '/clusters.txt', 'w') as f:
    f.write('# state ncells mass nselected cells\n'); [f.write('%d %d %.0f %d %s\n' % o) for o in out]
np.savez_compressed(a.out + '/cluster_targets.npz', T=np.array(T).astype(np.complex64), band=H['band'], n3=H['n3'], states=np.array([o[0] for o in out])); print('%d states -> %s' % (len(out), a.out))
