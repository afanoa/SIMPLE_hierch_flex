#!/usr/bin/env python3
"""Truth readout of one E-step iteration (problem 2): leaf labels and poses of a print_project_field dump against the simulation's
truth (gt_particles: conformer; all.txt: true Euler angles). NMI(conformer; leaf), conformer purity, swing eta2 and extent (scorecard
definitions); pose error = geodesic angle between refined and true rotations (SIMPLE euler2m: Rz(e3) Ry(e2) Rz(e1), rotmat signs),
median / mean / share above 10 deg over all non-junk particles and over those that sat in the NSMALL smallest leaves of the
reference labelling (--ref0); leaf changes against --prev. Writes the dump's labels to --out.
usage: p2_metrics.py DUMP --gt gt_particles.txt --true all.txt --ref0 labels_it1.txt [--prev labels.txt --out labels_next.txt --nsmall 10 --tag T]"""
import re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('dump'); ap.add_argument('--gt', required=True); ap.add_argument('--true', required=True); ap.add_argument('--ref0', required=True)
ap.add_argument('--prev', default=''); ap.add_argument('--out', default=''); ap.add_argument('--nsmall', type=int, default=10); ap.add_argument('--tag', default=''); a = ap.parse_args()
def oris(fn, keys):
    pat = {k: re.compile(r'(?:^| )%s=([-0-9.eE+]+)' % k) for k in keys}; rows = []
    for l in open(fn):
        if l.strip(): rows.append([float(pat[k].search(l).group(1)) if pat[k].search(l) else np.nan for k in keys])
    return np.array(rows)
X = oris(a.dump, ('state', 'e1', 'e2', 'e3', 'x', 'y')); Tr = oris(a.true, ('e1', 'e2', 'e3')); lab = X[:, 0].astype(int); nr = len(lab)
G = np.loadtxt(a.gt, comments='#'); conf = np.full(nr, -1); gi = G[:, 0].astype(int) - 1; v = (gi >= 0) & (gi < nr); conf[gi[v]] = G[v, 1].astype(int)
if a.out: np.savetxt(a.out, lab, fmt='%d')
def rz(t): c, s = np.cos(t), np.sin(t); z = np.zeros_like(t); o = np.ones_like(t); return np.stack([np.stack([c, s, z], -1), np.stack([-s, c, z], -1), np.stack([z, z, o], -1)], -2)
def ry(t): c, s = np.cos(t), np.sin(t); z = np.zeros_like(t); o = np.ones_like(t); return np.stack([np.stack([c, z, -s], -1), np.stack([z, o, z], -1), np.stack([s, z, c], -1)], -2)
def rot(e): e = np.deg2rad(e); return rz(e[:, 2]) @ ry(e[:, 1]) @ rz(e[:, 0])
Rr, Rt = rot(X[:, 1:4]), rot(Tr[:, :3]); ang = np.rad2deg(np.arccos(np.clip((np.einsum('nij,nij->n', Rr, Rt) - 1) / 2, -1, 1)))
on = (lab > 0) & (conf >= 0); L, C = lab[on], conf[on]
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); return I / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
_, Ci = np.unique(C, return_inverse=True); pur = sum(np.bincount(Ci[L == s]).max() for s in np.unique(L)) / len(L)
sw = (C >= 200) & (C <= 230); xs, ls = C[sw] - 200., L[sw]; mu = xs.mean(); eta = sum((ls == s).sum() * (xs[ls == s].mean() - mu) ** 2 for s in np.unique(ls)) / ((xs - mu) ** 2).sum()
ms = [xs[ls == s].mean() for s in np.unique(ls) if (ls == s).sum() >= 200 and (ls == s).sum() >= 0.5 * (L == s).sum()]
ref0 = np.loadtxt(a.ref0, dtype=int); cnt0 = {s: (ref0 == s).sum() for s in np.unique(ref0[ref0 > 0])}; small = sorted(cnt0, key=lambda s: cnt0[s])[:a.nsmall]; ins = np.isin(ref0, small)
good = conf > 0
def pst(m): q = ang[m & good]; return 'median %.2f mean %.2f >10deg %.1f %%' % (np.median(q), q.mean(), 100. * (q > 10).mean())
out = ['%s%d leaves over %d particles; NMI(conformer; leaf) %.3f, purity %.3f | swing eta2 %.3f, extent %s' % ((a.tag + ': ') if a.tag else '', len(np.unique(L)), len(L), nmi(C, L), pur, eta, ('%.1f .. %.1f in %d leaves' % (min(ms), max(ms), len(ms))) if ms else '-'),
       '   pose error vs truth (non-junk): all %s | %d smallest it1 leaves (%d particles) %s' % (pst(np.ones(nr, bool)), len(small), (ins & good).sum(), pst(ins))]
if a.prev:
    pv = np.loadtxt(a.prev, dtype=int); out.append('   leaf changes vs previous: %.1f %% of all, %.1f %% of the small leaves\' particles; small leaves now %s (it1 %s)'
               % (100. * (pv != lab).mean(), 100. * (pv[ins] != lab[ins]).mean(), ' '.join(str((lab == s).sum()) for s in small), ' '.join(str(cnt0[s]) for s in small)))
print('\n'.join(out))
