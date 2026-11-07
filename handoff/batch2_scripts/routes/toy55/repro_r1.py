#!/usr/bin/env python3
"""Run-to-run reproducibility of tree-node partitions from existing runs (no new refinement). Label files (one state per master row, 0 =
outside) for every run and step of node f692st r1 (script rflex13, its earlier cheap-mask run old_cheap20, the Fortran port) and of unit 7
(f692u7, f692u7n13f, f692u7b13f): within-run drift between consecutive steps, across-run agreement at the same step (NMI, ARI), group sizes,
and the closed tracer (old f692r closed states 15+16): per partition its most closed-enriched group and how that group's particles land in
the other runs.  usage: repro_r1.py LABELDIR > report.txt"""
import sys, os, itertools, numpy as np
L = sys.argv[1]; K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); oc = (old == 15) | (old == 16)
def lab(n): return np.loadtxt(os.path.join(L, n + '.txt'), dtype=int)
def nmi_ari(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((len(ua), len(ub))); np.add.at(C, (ia, ib), 1); n = C.sum(); P = C / n; pa, pb = P.sum(1), P.sum(0)
    nz = P > 0; I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); H = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
    c2 = lambda x: x * (x - 1) / 2.; sij, sa, sb = c2(C).sum(), c2(C.sum(1)).sum(), c2(C.sum(0)).sum(); e = sa * sb / c2(n)
    return I / np.sqrt(max(H(pa) * H(pb), 1e-12)), (sij - e) / (0.5 * (sa + sb) - e)
def top_closed(a, rows):
    g, enr = None, -1.
    base = oc[rows].mean()
    for s in np.unique(a[rows]):
        m = rows[a[rows] == s]
        if len(m) < 1000: continue
        e = oc[m].mean() / max(base, 1e-9)
        if e > enr: g, enr = s, e
    m = rows[a[rows] == g]; return g, m, enr, oc[m].sum() / max(oc[rows].sum(), 1)
SETS = {'f692st r1': ('S_B', [('script', ['S_flex', 'S_fflex', 'S_A', 'S_B']), ('script old_cheap20 (other mask)', [None, 'O_fflex', 'O_A', None]),
                              ('port', ['P_flex', 'P_fflex', 'P_A', 'P_B'])]),
        'unit 7': ('U_B', [('u7 (cheap mask)', ['U_flex', 'U_fflex', 'U_A', 'U_B']), ('u7n13f (falloff mask)', ['N_flex', 'N_fflex', 'N_A', 'N_B']),
                           ('u7b13f (same stage A as n13f, borrowed stage B)', [None, None, 'Bb_A', 'Bb_B'])])}
STEPS = ['whole flex', 'focused flex (stage-A seeds)', 'stage A', 'stage B']
for title, (rowsrc, runs) in SETS.items():
    rows = np.where(lab(rowsrc) > 0)[0]; print('=' * 100); print('%s: %d particles, old closed %d (%.1f %%)' % (title, len(rows), oc[rows].sum(), 100. * oc[rows].mean()))
    lb = {}
    for _, names in runs:
        for nm in names:
            if nm and os.path.exists(os.path.join(L, nm + '.txt')):
                a = lab(nm); lb[nm] = a; outside = int((a[rows] == 0).sum())
                sz = np.bincount(a[rows])[1:]; sz = sz[sz > 0]
                print('  %-7s %d groups, sizes %s%s' % (nm, len(sz), ' '.join(str(int(x)) for x in sorted(sz, reverse=True)), ('  (%d rows unlabelled)' % outside) if outside else ''))
    print('-- within-run drift between consecutive steps (NMI / ARI; 1 = unchanged)')
    for run, names in runs:
        out = []
        for i in range(3):
            a, b = names[i], names[i + 1]
            if a in lb and b in lb:
                m = rows[(lb[a][rows] > 0) & (lb[b][rows] > 0)]; v = nmi_ari(lb[a][m], lb[b][m]); out.append('%s -> %s %.2f / %.2f' % (STEPS[i], STEPS[i + 1], *v))
        if out: print('  %-48s %s' % (run, ';  '.join(out)))
    print('-- across runs at the same step (NMI / ARI)')
    for (r1, n1), (r2, n2) in itertools.combinations(runs, 2):
        out = []
        for i in range(4):
            a, b = n1[i], n2[i]
            if a in lb and b in lb:
                m = rows[(lb[a][rows] > 0) & (lb[b][rows] > 0)]; v = nmi_ari(lb[a][m], lb[b][m]); out.append('%s %.2f / %.2f' % (STEPS[i], *v))
        if out: print('  %s  vs  %s:  %s' % (r1, r2, ';  '.join(out)))
    print('-- closed tracer: most closed-enriched group per partition, and where its particles land in the other partitions of the same step')
    for i in range(4):
        parts = [(run, names[i]) for run, names in runs if names[i] in lb]
        if not parts: continue
        tops = {nm: top_closed(lb[nm], rows) for _, nm in parts}
        for run, nm in parts:
            g, m, enr, share = tops[nm]
            land = []
            for run2, nm2 in parts:
                if nm2 == nm: continue
                g2, m2, _, _ = tops[nm2]; c = np.bincount(lb[nm2][m], minlength=lb[nm2].max() + 1)[1:]
                land.append('in %s: %.0f %% in its closed group, largest share %.0f %%, Jaccard %.2f' % (nm2, 100. * c[g2 - 1] / len(m), 100. * c.max() / len(m),
                            len(np.intersect1d(m, m2)) / len(np.union1d(m, m2))))
            print('  [%s] %-7s group %d: %d particles, %.1fx enriched, %.0f %% of the node\'s old closed%s' % (STEPS[i], nm, g, len(m), enr, 100. * share,
                  ('; ' + '; '.join(land)) if land else ''))
