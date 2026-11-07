#!/usr/bin/env python3
"""usage: validate_tree.py NAME TREE_DIR TRUTH KIND [--tmasks A.mrc,B.mrc,C.mrc] [--table factorial/table.txt (m8c19 on; older trees: factorial_v3/table.txt)]
Validation of a port hetero_tree (TREE_DIR = .../1_hetero_tree) against toy truth, one summary line (prefixed SUMMARY) + details.
Both kinds also report the junk kept in leaves (overall, worst leaf). KIND sites (simsites truth: junk, A_domain, B_blob, C_closed,
combo): junk precision / recall of the removed rows (label 0); leaf purity (genuine; and with junk as its own class) =
size-weighted share of each leaf's dominant combination (genuine particles); A split = share of genuine particles in A-pure leaves
(A-present share <= 0.2 or >= 0.8) and the A-mixed leaves; the tree's most closed leaf (precision / recall); the per-site layer's best
closed pole (labels: precision / recall; posterior AUC) from TABLE.
KIND arms (IgM truth: junk + arm angles): junk precision / recall; eta^2 of each angle over the leaves; per angle the best site posterior
|Spearman| from TABLE."""
import os, sys, numpy as np
NAME, D, TF, KIND = sys.argv[1:5]; opt = dict(zip(sys.argv[5::2], sys.argv[6::2])); TAB = opt.get('--table', 'factorial/table.txt')
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L)
def col(c): v = np.full(NM, np.nan); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; rem = L == 0; g = (L > 0) & ~junk; leaves = [k for k in np.unique(L) if k > 0]
names = [l.split()[0] for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]
jp, jr = (rem & junk).sum() / max(rem.sum(), 1), (rem & junk).sum() / max(junk.sum(), 1)
inl = L > 0; jl = {k: junk[L == k].mean() for k in leaves}; jin = (junk & inl).sum() / max(inl.sum(), 1)   # junk kept in leaves (empty-only rule)
JSUM = 'junk in leaves %.3f (max leaf %.3f, %d leaves > 0.2)' % (jin, max(jl.values()) if jl else 0., sum(v > 0.2 for v in jl.values()))
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def spr(x, y): return float(np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1])
def table_rows():
    p = TAB if os.path.isabs(TAB) else os.path.join(D, TAB)
    if not os.path.exists(p): return []
    out = []
    for w in [l.split() for l in open(p) if l.strip() and not l.startswith('#')]:
        Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T; r = np.loadtxt(w[4], dtype=int).ravel()
        rows = (r > 0) if (len(r) == NM and r.max() <= 1) else np.isin(np.arange(NM), r - 1); out.append((w[0], np.loadtxt(w[2], dtype=int), Q, rows))
    return out
print('== %s: %s, %d leaves, %d removed; junk removal precision %.3f recall %.3f' % (NAME, D, len(leaves), rem.sum(), jp, jr))
if KIND == 'sites':
    A, C, cmb = col('A_domain') > 0, col('C_closed') > 0, col('combo')
    pur = sum(max(np.bincount(cmb[(L == k) & g].astype(int), minlength=9)) for k in leaves) / max(g.sum(), 1)
    cj = np.where(junk, 0, np.nan_to_num(cmb)).astype(int); pura = sum(max(np.bincount(cj[L == k], minlength=9)) for k in leaves) / max(inl.sum(), 1)
    ash = {k: A[(L == k) & g].mean() for k in leaves}; apure = sum(((L == k) & g).sum() for k in leaves if ash[k] <= 0.2 or ash[k] >= 0.8) / max(g.sum(), 1)
    cl = {k: ((C & (L == k) & g).sum() / max(((L == k) & g).sum(), 1), (C & (L == k) & g).sum() / max((C & g).sum(), 1), ((L == k) & g).sum()) for k in leaves}
    kb = max(cl, key=lambda k: cl[k][0])
    for k in leaves: print('   leaf %2d %-14s %7d genuine, junk %.3f | domain A %.2f%s | closed %.3f (recall %.3f)' % (k, names[k - 1] if k <= len(names) else '?', cl[k][2], jl[k], ash[k], '' if ash[k] <= 0.2 or ash[k] >= 0.8 else '  <- A-MIXED', cl[k][0], cl[k][1]))
    best = None
    for nm, lab, Q, rows in table_rows():
        gg = rows & ~junk
        for j in range(Q.shape[1]):
            m = (lab == j + 1) & rows; P = Q[:, j]; ok = gg & (P >= 0)
            if not m.any() or C[ok].sum() == 0: continue
            pr, rc = (m & C).sum() / max((m & ~junk).sum(), 1), (m & C).sum() / max((gg & C).sum(), 1); a_ = auc(P[ok], C[ok])
            if best is None or pr > best[0]: best = (pr, rc, max(a_, 1 - a_), '%s pole %d' % (nm, j + 1), m.sum())
    site = 'site layer best closed pole %s: %d labelled, precision %.3f recall %.3f, posterior AUC %.3f' % (best[3], best[4], best[0], best[1], best[2]) if best else 'site layer: no table'
    print('SUMMARY %-12s junk %.3f/%.3f | %s | %d leaves, purity %.3f (incl. junk %.3f) | A-pure share %.3f (%d A-mixed leaves) | tree closed leaf %.3f/%.3f | %s'
          % (NAME, jp, jr, JSUM, len(leaves), pur, pura, apure, sum(1 for k in leaves if 0.2 < ash[k] < 0.8), cl[kb][0], cl[kb][1], site))
else:
    cols = [c for c in hd[3:]]
    def eta2(lab, y): mu = y.mean(); return sum(((lab == k).sum() * (y[lab == k].mean() - mu) ** 2) for k in np.unique(lab)) / max(((y - mu) ** 2).sum(), 1e-30)
    e = {c: eta2(L[g], col(c)[g]) for c in cols}; bs = {}
    for nm, lab, Q, rows in table_rows():
        gg = rows & ~junk
        for c in cols:
            y = col(c)
            for j in range(Q.shape[1]):
                P = Q[:, j]; ok = gg & (P >= 0)
                if ok.sum() < 100: continue
                r_ = abs(spr(P[ok], y[ok]))
                if r_ > bs.get(c, (0., ''))[0]: bs[c] = (r_, nm)
    print('SUMMARY %-12s junk %.3f/%.3f | %s | %d leaves | leaf eta^2 %s | site |rho| %s' % (NAME, jp, jr, JSUM, len(leaves), ' '.join('%s %.2f' % (c, e[c]) for c in cols),
          ' '.join('%s %.2f' % (c, bs[c][0]) if c in bs else '%s -' % c for c in cols)))
