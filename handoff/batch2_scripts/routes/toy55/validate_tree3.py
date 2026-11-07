#!/usr/bin/env python3
"""usage: validate_tree3.py NAME TREE_DIR TRUTH KIND [--table factorial/table.txt]
validate_tree2.py plus junk isolation, NMI and the truth support of every split.
- Removed = rows in no leaf plus the leaves the tree names 'junk'; junk-dominated leaves (junk share > 0.5) are junk leaves; junk in
  genuine leaves is counted over the rest.
- Junk isolation per junk group (sites: jtype; arms: one group): R = share of the group removed or in junk leaves; P = junk share of
  the junk buckets (the removed rows, each junk leaf) where the group is the larger junk group.
- NMI(partition; truth class): the leaves plus the removed rows as one cell; class = combo, each junk group its own class.
- Splits: every inner node (children = the next node or state leaf on each leaf's path) with its children's sizes and the truth
  factor that differs most between children of >= 200 rows (sites: share difference of junk, each junk group, A, B, C; arms: junk share difference,
  eta^2 of the arm angles): supported >= 0.2 (share) or 0.1 (eta^2), weak >= 0.1 / 0.03, else none (a null-type split)."""
import os, sys, numpy as np
NAME, D, TF, KIND = sys.argv[1:5]; opt = dict(zip(sys.argv[5::2], sys.argv[6::2])); TAB = opt.get('--table', 'factorial/table.txt')
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1
L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); NM = len(L)
def col(c): v = np.full(NM, np.nan); v[tr] = T[:, hd.index(c)]; return v
LV = [l.split() for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]; names = [w[0] for w in LV]
junk = col('junk') > 0; jt = np.nan_to_num(col('jtype')).astype(int) if 'jtype' in hd else junk.astype(int)
leaves = [k for k in np.unique(L) if k > 0]; jmark = [k for k in leaves if k <= len(LV) and LV[k - 1][1] == 'junk']
rem = (L == 0) | np.isin(L, jmark); glv = [k for k in leaves if k not in jmark]
jdom = [k for k in glv if junk[L == k].mean() > 0.5]; gen = [k for k in glv if k not in jdom]
g = ~junk & np.isin(L, glv); iso = rem | np.isin(L, jdom)
jp, jr = (rem & junk).sum() / max(rem.sum(), 1), (rem & junk).sum() / max(junk.sum(), 1)
jl = {k: junk[L == k].mean() for k in leaves}; ing = np.isin(L, gen); jin = (junk & ing).sum() / max(ing.sum(), 1)
JSUM = 'junk in genuine leaves %.3f (worst %.3f, %d > 0.2; %d junk leaves)' % (jin, max([jl[k] for k in gen], default=0.), sum(jl[k] > 0.2 for k in gen), len(jmark) + len(jdom))
bk = [L == 0] + [L == k for k in jmark + jdom]; JG = []; qs = sorted(np.unique(jt[junk]))
def jcount(b, q): return ((jt == q) & b & junk).sum()
for q in qs:
    own = [b for b in bk if jcount(b, q) > 0 and jcount(b, q) == max(jcount(b, r) for r in qs)]
    bm = np.any(own, axis=0) if own else np.zeros(NM, bool)
    JG.append('j%d %s/%.3f' % (q, '%.3f' % (junk[bm].mean()) if bm.any() else 'n/a', (iso & (jt == q) & junk).sum() / max(((jt == q) & junk).sum(), 1)))
def nmi(a, b):
    _, ia = np.unique(a, return_inverse=True); _, ib = np.unique(b, return_inverse=True); n = len(a)
    P = np.zeros((ia.max() + 1, ib.max() + 1)); np.add.at(P, (ia, ib), 1. / n); pa, pb = P.sum(1), P.sum(0); nz = P > 0
    I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); H = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
    return 2 * I / max(H(pa) + H(pb), 1e-12)
cls = np.where(junk, -jt, np.nan_to_num(col('combo'))).astype(int); NMI = nmi(L[tr], cls[tr])
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def spr(x, y): return float(np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1])
def eta2(lab, y): mu = y.mean(); return sum(((lab == k).sum() * (y[lab == k].mean() - mu) ** 2) for k in np.unique(lab)) / max(((y - mu) ** 2).sum(), 1e-30)
def table_rows():
    p = TAB if os.path.isabs(TAB) else os.path.join(D, TAB)
    if not os.path.exists(p): return []
    out = []
    for w in [l.split() for l in open(p) if l.strip() and not l.startswith('#')]:
        Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T; r = np.loadtxt(w[4], dtype=int).ravel()
        rows = (r > 0) if (len(r) == NM and r.max() <= 1) else np.isin(np.arange(NM), r - 1); out.append((w[0], np.loadtxt(w[2], dtype=int), Q, rows))
    return out
def path(nm):
    node, _, rest = nm.partition('_'); d = node[1:]; p = ['r' + d[:i] for i in range(len(d) + 1)]
    return p + [nm] if rest.startswith('s') else p
inner = {}
for k in leaves:
    p = path(names[k - 1]) if k <= len(names) else ['r', '?%d' % k]
    for i, x in enumerate(p[:-1]): inner.setdefault(x, {}).setdefault(p[i + 1], []).append(k)
ang = [c for c in hd[3:]] if KIND == 'arms' else []; ht = ~np.isnan(col('junk'))
SPL, NS = [], {'supported': 0, 'weak': 0, 'none': 0}
for x in sorted(inner, key=lambda s: (len(s), s)):
    ch = inner[x]
    if len(ch) < 2: continue
    cm = {c: np.isin(L, ks) for c, ks in ch.items()}; big = [c for c in ch if cm[c].sum() >= 200]
    eff = {'junk': np.ptp([junk[cm[c]].mean() for c in big]) if len(big) > 1 else 0.}
    for q in (qs if len(qs) > 1 else []): eff['j%d' % q] = np.ptp([(junk & (jt == q))[cm[c]].mean() for c in big]) if len(big) > 1 else 0.
    if KIND == 'sites':
        for f in ('A_domain', 'B_blob', 'C_closed'):
            gb = [c for c in big if (cm[c] & g).sum() >= 200]; v = col(f) > 0
            eff[f[0]] = np.ptp([v[cm[c] & g].mean() for c in gb]) if len(gb) > 1 else 0.
        top = max(eff, key=eff.get); val = eff[top]; sup, wk = val >= 0.2, val >= 0.1
    else:
        nx = np.isin(L, sum(ch.values(), [])) & g; lab = np.zeros(NM, int)
        for i, c in enumerate(ch): lab[cm[c]] = i
        for c in ang: eff[c] = eta2(lab[nx], col(c)[nx]) if nx.sum() > 100 else 0.
        ea = max(ang, key=lambda c: eff[c]); sup, wk = eff['junk'] >= 0.2 or eff[ea] >= 0.1, eff['junk'] >= 0.1 or eff[ea] >= 0.03
        top = 'junk' if eff['junk'] / 0.2 >= eff[ea] / 0.1 else ea; val = eff[top]
    st = 'supported' if sup else 'weak' if wk else 'none'; NS[st] += 1; nxa = np.isin(L, sum(ch.values(), [])) & ht
    SPL.append('%s K%d %s -> %s (%s %.2f; NMI %.3f)' % (x, len(ch), '/'.join('%d' % cm[c].sum() for c in ch), st, top, val, nmi(L[nxa], cls[nxa])))
print('== %s: %s, %d leaves (%d named junk), %d removed; junk removal precision %.3f recall %.3f; NMI %.3f' % (NAME, D, len(leaves), len(jmark), rem.sum(), jp, jr, NMI))
for s in SPL: print('   split ' + s)
SSUM = 'splits %d supported, %d weak, %d none%s' % (NS['supported'], NS['weak'], NS['none'], ''.join(' [%s]' % s.split(' (')[0] for s in SPL if '-> none' in s or '-> weak' in s))
if KIND == 'sites':
    A, C, cmb = col('A_domain') > 0, col('C_closed') > 0, col('combo')
    pur = sum(max(np.bincount(cmb[(L == k) & g].astype(int), minlength=9)) for k in glv) / max(g.sum(), 1)
    cj = np.where(junk, 0, np.nan_to_num(cmb)).astype(int); inl = L > 0; pura = sum(max(np.bincount(cj[L == k], minlength=9)) for k in leaves) / max(inl.sum(), 1)
    ash = {k: A[(L == k) & g].mean() if ((L == k) & g).any() else 0. for k in glv}; apure = sum(((L == k) & g).sum() for k in glv if ash[k] <= 0.2 or ash[k] >= 0.8) / max(g.sum(), 1)
    cl = {k: ((C & (L == k) & g).sum() / max(((L == k) & g).sum(), 1), (C & (L == k) & g).sum() / max((C & g).sum(), 1), ((L == k) & g).sum()) for k in glv}
    kb = max(cl, key=lambda k: cl[k][0])
    for k in leaves: print('   leaf %2d %-14s %7d rows, junk %.3f%s | domain A %.2f%s | closed %.3f (recall %.3f)' % (k, names[k - 1] if k <= len(names) else '?', (L == k).sum(), jl[k], ' NAMED JUNK' if k in jmark else ' junk leaf' if k in jdom else '', ash.get(k, 0.), '  <- A-MIXED' if k in gen and 0.2 < ash[k] < 0.8 else '', cl[k][0] if k in cl else 0., cl[k][1] if k in cl else 0.))
    best = None
    for nm, lab, Q, rows in table_rows():
        gg = rows & ~junk
        for j in range(Q.shape[1]):
            m = (lab == j + 1) & rows; P = Q[:, j]; ok = gg & (P >= 0)
            if not m.any() or C[ok].sum() == 0: continue
            pr, rc = (m & C).sum() / max((m & ~junk).sum(), 1), (m & C).sum() / max((gg & C).sum(), 1); a_ = auc(P[ok], C[ok])
            if best is None or pr > best[0]: best = (pr, rc, max(a_, 1 - a_), '%s pole %d' % (nm, j + 1), m.sum())
    site = 'site layer best closed pole %s: %d labelled, precision %.3f recall %.3f, posterior AUC %.3f' % (best[3], best[4], best[0], best[1], best[2]) if best else 'site layer: no table'
    print('SUMMARY %-12s junk %.3f/%.3f | junk groups P/R %s | %s | %d leaves, purity %.3f (incl. junk %.3f), NMI %.3f | A-pure share %.3f (%d A-mixed leaves) | tree closed leaf %.3f/%.3f | %s | %s'
          % (NAME, jp, jr, ' '.join(JG), JSUM, len(leaves), pur, pura, NMI, apure, sum(1 for k in gen if 0.2 < ash[k] < 0.8), cl[kb][0], cl[kb][1], site, SSUM))
else:
    e = {c: eta2(L[g], col(c)[g]) for c in ang}; bs = {}
    for k in leaves: print('   leaf %2d %-14s %7d rows, junk %.3f%s' % (k, names[k - 1] if k <= len(names) else '?', (L == k).sum(), jl[k], ' NAMED JUNK' if k in jmark else ' junk leaf' if k in jdom else ''))
    for nm, lab, Q, rows in table_rows():
        gg = rows & ~junk
        for c in ang:
            y = col(c)
            for j in range(Q.shape[1]):
                P = Q[:, j]; ok = gg & (P >= 0)
                if ok.sum() < 100: continue
                r_ = abs(spr(P[ok], y[ok]))
                if r_ > bs.get(c, (0., ''))[0]: bs[c] = (r_, nm)
    print('SUMMARY %-12s junk %.3f/%.3f | junk groups P/R %s | %s | %d leaves, NMI %.3f | leaf eta^2 %s | site |rho| %s | %s' % (NAME, jp, jr, ' '.join(JG), JSUM, len(leaves), NMI,
          ' '.join('%s %.2f' % (c, e[c]) for c in ang), ' '.join('%s %.2f' % (c, bs[c][0]) if c in bs else '%s -' % c for c in ang), SSUM))
