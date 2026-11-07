#!/usr/bin/env python3
"""Is a run's organization meaningful? Truth readout on a synthetic set of composition classes (apo conformers < 100, holo
101-199, junk 0) and CONTINUOUS families (FAMS: name:first_id:A_per_step, e.g. swing:200:1.0 = conformer 200+k is step k,
k A along the motion). Per continuous family:
(1) top split: each family of the run (expand/family<F>_subset.txt) by class, NMI(class; family);
(2) branch coherence: share of the family's particles (all, first and last QEND of its steps) in their top run-family,
    top region board and top leaf;
(3) leaf purity: leaves (level-1 states of region / zoom boards; first board holding a particle) with >= 80 % of one
    continuous family, and the spread of their steps (A) against the family's whole spread;
(4) coverage: the motion range spanned by those pure leaves (min / max mean position), against the true range;
(5) the same for every consolidation of the run (states.txt): final states.
usage: org_readout.py RUN GT [--fams swing:200:1.0,breath:300:0.7 --qend 0.1]"""
import sys, os, glob, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('gt'); ap.add_argument('--fams', default='swing:200:1.0'); ap.add_argument('--qend', type=float, default=0.1); a = ap.parse_args()
G = np.loadtxt(a.gt, comments='#'); conf = G[:, 1].astype(int); n = len(conf)
FAMS = [(w.split(':')[0], int(w.split(':')[1]), float(w.split(':')[2])) for w in a.fams.split(',')]; FAMS.sort(key=lambda f: f[1])
cls = np.where(conf == 0, 0, np.where(conf >= 101, 2, 1)); names = {0: 'junk', 1: 'apo', 2: 'holo'}; pos = {}
for i, (nm, i0, st) in enumerate(FAMS):
    hi = FAMS[i + 1][1] if i + 1 < len(FAMS) else 10 ** 9; m = (conf >= i0) & (conf < hi); cls[m] = 3 + i; names[3 + i] = nm; pos[3 + i] = np.where(m, (conf - i0) * st, np.nan)
def nmi(x, y):
    xi, yi = np.unique(x, return_inverse=True)[1], np.unique(y, return_inverse=True)[1]; J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); P = J / J.sum()
    px, py = P.sum(1), P.sum(0); nz = P > 0; I = (P[nz] * np.log(P[nz] / np.outer(px, py)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
    return I / max(np.sqrt(h(px) * h(py)), 1e-12)
def top_share(lab, m): v = lab[m]; v = v[v > 0]; return (np.bincount(v).max() / m.sum()) if len(v) and m.sum() else 0.
def ends(c): p = pos[c]; t = np.nanmax(p); return (cls == c) & (p <= a.qend * t), (cls == c) & (p >= (1 - a.qend) * t), t
fam = np.zeros(n, int)
for F in range(1, 9):
    f = '%s/expand/family%d_subset.txt' % (a.run, F)
    if os.path.isfile(f): fam[np.loadtxt(f, dtype=int).ravel() - 1] = F
print('== (1) top split (run families)')
for F in sorted(set(fam) - {0}):
    m = fam == F; print('  family %d: %d particles = %s' % (F, m.sum(), ', '.join('%s %.0f %%' % (names[c], 100 * np.mean(cls[m] == c)) for c in sorted(names))))
inf = fam > 0; comp = np.where(cls >= 3, 1, cls)
print('  NMI(class; family) %.3f   NMI(composition apo/holo/junk; family) %.3f' % (nmi(cls[inf], fam[inf]), nmi(comp[inf], fam[inf])))
def leaves(level):
    L = {}
    for f in sorted(glob.glob('%s/fam?/%s*/level1/states_clusters.txt' % (a.run, level))):
        s = np.loadtxt(f, dtype=int)
        if len(s) == n: L[f.split('/')[-4] + '/' + f.split('/')[-3]] = s
    board = np.zeros(n, int); leaf = np.zeros(n, int); nl = 0
    for i, (b, s) in enumerate(L.items()):   # a particle's board / leaf = the first one holding it (boards of a family overlap)
        new = (s > 0) & (leaf == 0)
        for k in np.unique(s[new]): nl += 1; leaf[new & (s == k)] = nl
        board[new] = i + 1
    return len(L), board, leaf
def purity_report(tag, lab):
    for c in sorted(pos):
        lo, hi, tmax = ends(c); mc = cls == c; rows = []
        for k in np.unique(lab[lab > 0]):
            m = lab == k
            if m.sum() < 200: continue
            sh = np.mean(cls[m] == c); p = pos[c][m & mc]; rows.append((sh, m.sum(), np.nanmean(p) if len(p) else np.nan, np.nanstd(p) if len(p) > 1 else np.nan))
        pure = sorted([r for r in rows if r[0] >= 0.8], key=lambda r: r[2])
        print('  %s %s: top %.2f (ends %.2f / %.2f); %d pure states%s' % (tag, names[c], top_share(lab, mc), top_share(lab, lo), top_share(lab, hi), len(pure),
              ('; span %.1f .. %.1f of 0 .. %.0f A, spread median %.1f A (whole %.1f A): %s' % (pure[0][2], pure[-1][2], tmax, np.median([r[3] for r in pure]), np.nanstd(pos[c][mc]), ' '.join('%.0f+-%.0f(%d)' % (r[2], r[3], r[1]) for r in pure))) if pure else ''))
print('== (2-4) continuous families: branch coherence, pure leaves, coverage')
for c in sorted(pos):
    lo, hi, _ = ends(c); print('  %s: in its top run-family %.2f (ends %.2f / %.2f)' % (names[c], top_share(fam, cls == c), top_share(fam, lo), top_share(fam, hi)))
for lv in ('region', 'zoomregion'):
    nb, board, leaf = leaves(lv)
    if not nb: continue
    for c in sorted(pos):
        lo, hi, _ = ends(c); print('  %s level (%d boards) %s: held %.2f, top board %.2f (ends %.2f / %.2f)' % (lv, nb, names[c], np.mean(leaf[cls == c] > 0), top_share(board, cls == c), top_share(board, lo), top_share(board, hi)))
    purity_report('%s leaves' % lv, leaf)
print('== (5) consolidated states')
for C in sorted(glob.glob(a.run + '/consolidated_*')) + sorted(glob.glob(a.run + '/mergeup/consolidated*')) + sorted(glob.glob(a.run + '/landscape*')) + sorted(glob.glob(a.run + '/graph_*')):
    f = C + '/states.txt'
    if os.path.isfile(f):
        st = np.loadtxt(f, dtype=int)
        if len(st) == n: print(' %s (%d states)' % (os.path.relpath(C, a.run), len(np.unique(st[st > 0])))); purity_report('  ', st)
