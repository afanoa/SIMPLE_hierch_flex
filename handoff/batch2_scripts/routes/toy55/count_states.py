#!/usr/bin/env python3
"""How many DISTINCT, REPRODUCIBLE states does a set of klem boards resolve? Cells survive by the cross-half
reproducibility of their deviation (klem_states.txt, >= --minrep); surviving cells from all boards are clustered by the
cross-half similarity of their deviations (corr of one half's deviation of cell i with the other half's of cell j,
symmetrised, from klem_halves.npz, <= 12 A), single linkage at --sim thresholds; a cluster = one state. Optional truth
(--gt) gives the dominant conformer per cluster and the number of conformers represented; optional deposited matches
(per-board cell_matches.txt) give the deposited conformer per cluster. usage: count_states.py DIR1 DIR2 ... [--minrep 0.5] [--sim 0.7,0.8,0.9] [--gt gt.txt]"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('dirs', nargs='+'); ap.add_argument('--minrep', type=float, default=0.5); ap.add_argument('--sim', default='0.6,0.7,0.8,0.9'); ap.add_argument('--gt', default=''); ap.add_argument('--save', default='', help='write per-cell cluster labels at --at similarity: board_index dir cell cluster weight dominant_deposited corr'); ap.add_argument('--at', type=float, default=0.6); a = ap.parse_args()
D0 = []; D1 = []; W = []; src = []; cell = []; rep = []; conf_of_cell = []; dep_of_cell = []
gt = {}
if a.gt:
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): t = l.split(); gt[int(t[0])] = int(t[1])
for j, d in enumerate(a.dirs):
    if not os.path.exists(d + '/klem_halves.npz'): print('  %s: no klem_halves.npz (older run), skipped' % d); continue
    H = np.load(d + '/klem_halves.npz'); st = np.loadtxt(d + '/klem_states.txt', comments='#'); crep = st[:, 2]; ngrid = len(H['w'])
    L = np.loadtxt(d + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]
    dep = {}
    for f in (d + '/cell_matches.txt', d + '/matched_rec/cell_matches.txt', d + '/../matched_rec/cell_matches.txt', d + '/../cell_matches.txt'):
        if os.path.exists(f):
            for r in np.loadtxt(f, comments='#'): dep[int(r[0]) - 1] = (int(r[1]), float(r[2]))
            break
    for k in range(ngrid):
        if crep[k] < a.minrep: continue
        D0.append(H['D0'][k]); D1.append(H['D1'][k]); W.append(float(H['w'][k])); src.append(j); cell.append(k); rep.append(float(crep[k]))
        m = (gm > 0.5) & (am == k); c = np.array([gt.get(int(p), 0) for p in pind[m]]) if gt else np.zeros(0, int)
        conf_of_cell.append((int(np.bincount(c[c > 0]).argmax()), float(np.bincount(c[c > 0]).max() / max((c > 0).sum(), 1))) if (c > 0).any() else (0, 0.)); dep_of_cell.append(dep.get(k, (0, 0.)))
D0 = np.array(D0); D1 = np.array(D1); n = len(D0); print('%d surviving cells (reproducibility >= %.2f) from %d boards' % (n, a.minrep, len(a.dirs)))
if n == 0: raise SystemExit
nrm0 = np.sqrt(np.sum(np.abs(D0) ** 2, 1)); nrm1 = np.sqrt(np.sum(np.abs(D1) ** 2, 1)); C = np.real(D0 @ np.conj(D1).T) / np.maximum(np.outer(nrm0, nrm1), 1e-30); S = 0.5 * (C + C.T)
print('cross-half similarity between surviving cells: median %.2f, 90th pct %.2f; self-reproducibility median %.2f' % (np.median(S[np.triu_indices(n, 1)]), np.percentile(S[np.triu_indices(n, 1)], 90), np.median(np.diag(C))))
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
Dm = np.clip(1 - S, 0, 2); np.fill_diagonal(Dm, 0); Z = linkage(squareform(Dm, checks=False), method='average')
for s in [float(x) for x in a.sim.split(',')]:
    lab = fcluster(Z, t=1 - s, criterion='distance'); nk = lab.max(); sizes = np.array([np.sum(np.array(W)[lab == k]) for k in range(1, nk + 1)])
    line = '  similarity %.1f: %d states; particles per state median %d (min %d, max %d)' % (s, nk, np.median(sizes), sizes.min(), sizes.max())
    if gt:
        dom = [max(set(conf_of_cell[i][0] for i in np.where(lab == k)[0] if conf_of_cell[i][0] > 0), key=[conf_of_cell[i][0] for i in np.where(lab == k)[0]].count, default=0) for k in range(1, nk + 1)]
        line += ' | conformers represented %d/%d (states per conformer: %s)' % (len(set(x for x in dom if x > 0)), len(set(v for v in gt.values() if v > 0)), ' '.join('%d:%d' % (c, dom.count(c)) for c in sorted(set(x for x in dom if x > 0))))
    if any(v[0] > 0 for v in dep_of_cell):
        dd = [max(set(dep_of_cell[i][0] for i in np.where(lab == k)[0] if dep_of_cell[i][1] >= 0.3), key=[dep_of_cell[i][0] for i in np.where(lab == k)[0]].count, default=0) for k in range(1, nk + 1)]
        line += ' | deposited conformers identified %d (states per conformer: %s; unmatched states %d)' % (len(set(x for x in dd if x > 0)), ' '.join('%d:%d' % (c, dd.count(c)) for c in sorted(set(x for x in dd if x > 0))), dd.count(0))
    print(line)
    if a.save and abs(s - a.at) < 1e-9:
        with open(a.save, 'w') as f:
            f.write('# board dir cell cluster weight deposited_of_cell corr_of_cell\n')
            for i in range(n): f.write('%d %s %d %d %.1f %d %.3f\n' % (src[i], a.dirs[src[i]], cell[i], lab[i], W[i], dep_of_cell[i][0], dep_of_cell[i][1]))
        print('  cluster labels at %.1f -> %s' % (s, a.save))
