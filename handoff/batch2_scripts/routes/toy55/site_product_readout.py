#!/usr/bin/env python3
"""Offline site-product readout of a finished hetero_tree (no new compute): the kept sites pooled at the root (factorial/table.txt, level r)
are the candidate motions. GATE (the site layer's own numbers, factorial/sites.txt): a site is a motion only if its pole has >= MINPOP
hard members (posterior >= 0.5) and decisiveness = hard members / soft pole size >= 0.5 (a quality or junk tail, a duplicate and a null
site have soft, spread posteriors: 0.00-0.40 on every control; real and continuum sites 0.81-3.3). Per pair of gated sites the pole-label NMI (posterior >= 0.5) and the posterior correlation r; sites with NMI > 0.05 or
|r| > 0.2 are one motion (single linkage) -> the factor count; per particle the product of the motions' pole bits (the strongest site per
motion = largest pole-rule margin, factorial/sites.txt) -> product labels. Truth / field columns (optional) only score: per pole the best
category F1 of each categorical column (<= 20 integer values) and Spearman(posterior, column) of each continuous one.
usage: site_product_readout.py RUN_DIR(holding 1_hetero_tree) tag=.. [minpop=from run_manifest.txt keys, else 900] [truth=FILE] [cols=a,b,..] [fields=PRINT_PROJECT_FIELD_DUMP]
       [labels_out=FILE] (product label per project row)"""
import sys, os, re, itertools, numpy as np
from scipy.stats import spearmanr
args = dict(a.split('=', 1) if '=' in a else ('run', a) for a in sys.argv[1:]); R = args['run'].rstrip('/') + '/1_hetero_tree'; tag = args.get('tag', R)
def loc(p): return p if os.path.exists(p) else os.path.join(R, p[p.index('1_hetero_tree') + len('1_hetero_tree') + 1:])
lab = np.loadtxt(R + '/hetero_tree_labels.txt', dtype=int, ndmin=1); N = len(lab)
sites, margin, soft, hard = [], {}, {}, {}
manf = args['run'].rstrip('/') + '/run_manifest.txt'
mp = re.search(r'minpop=(\d+)', open(manf).read()) if os.path.exists(manf) else None
minpop = int(args.get('minpop', mp.group(1) if mp else 900))
if os.path.exists(R + '/factorial/sites.txt'):
    for l in open(R + '/factorial/sites.txt'):
        if l.startswith('#') or '|' not in l: continue
        w = l.split('|'); nm = w[0].split()[0]
        fz = lambda v: 1e4 if '*' in v else float(v)   # ******* = an overflowed print of a huge value
        try: margin[nm] = max(fz(x) for x in w[3].split()); soft[nm] = fz(w[2].split()[0]); hard[nm] = fz(w[4].split()[0])
        except (ValueError, IndexError): margin[nm] = 0.
for l in open(R + '/factorial/table.txt'):
    if l.startswith('#') or not l.strip(): continue
    nm, lev, lf, pf = l.split()[:4]
    if lev != 'r': continue
    pp = np.loadtxt(loc(pf), ndmin=2)
    if pp.shape[0] != N and pp.shape[1] == N: pp = pp.T
    sites.append((nm, pp[:, 0]))
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; p = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(p, (a, b), 1); p /= p.sum()
    pa, pb = p.sum(1), p.sum(0); nz = p > 0; mi = (p[nz] * np.log(p[nz] / np.outer(pa, pb)[nz])).sum()
    ha, hb = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(), -(pb[pb > 0] * np.log(pb[pb > 0])).sum(); return mi / max(np.sqrt(ha * hb), 1e-12)
# truth / fields
cols = {}
if 'truth' in args:
    hdr = [l for l in open(args['truth']) if l.startswith('#')][0].lstrip('#').split(); t = np.loadtxt(args['truth'], comments='#', ndmin=2)
    rows = t[:, hdr.index('row')].astype(int) if 'row' in hdr else np.arange(1, len(t) + 1)
    want = args.get('cols', ','.join(h for h in hdr if h != 'row')).split(',')
    for c in want:
        if c in hdr: v = np.full(N, np.nan); v[rows - 1] = t[:, hdr.index(c)]; cols[c] = v
if 'fields' in args:
    F = {k: [] for k in ('e1', 'e2', 'e3', 'dfx', 'x', 'y')}
    for l in open(args['fields']):
        d = dict(x.split('=', 1) for x in l.split() if '=' in x)
        for k in F: F[k].append(float(d.get(k, 'nan')))
    for k, v in F.items():
        if len(v) == N: cols['in_' + k] = np.array(v)
    if 'in_e1' in cols:   # view direction as a unit vector, so a pole can be scored against views
        th, ph = np.radians(cols['in_e2']), np.radians(cols['in_e1']); cols['in_vx'], cols['in_vy'], cols['in_vz'] = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
gen = np.ones(N, bool) if 'junk' not in cols else (cols['junk'] == 0)
print('== %s: %d rows, %d leaves, %d kept root-level site(s)' % (tag, N, lab.max(), len(sites)))
if not sites:
    print('FACTORS %s 0 (no root-level site: the product step adds nothing)' % tag); sys.exit(0)
for nm, p in sites:
    pole = p >= 0.5; out = []
    for c, v in cols.items():
        ok = np.isfinite(v) & (p >= 0)
        u = np.unique(v[ok])
        if len(u) <= 20 and np.allclose(u, np.round(u)):
            best = max((2. * (pole[ok] & (v[ok] == x)).sum() / (pole[ok].sum() + (v[ok] == x).sum()), x) for x in u)
            out.append('%s F1 %.2f (=%g)' % (c, best[0], best[1]))
        else:
            out.append('%s rho %+.2f' % (c, spearmanr(p[ok & gen], v[ok & gen])[0]))
    print('site %-8s pole %6d (%.3f) margin %7.1f | %s' % (nm, pole.sum(), pole.mean(), margin.get(nm, np.nan), '; '.join(out)))
keep = []
for nm, p in sites:
    nh = hard.get(nm, (p >= 0.5).mean()) * N; dec = nh / soft[nm] if soft.get(nm) else 0.
    ok = nh >= minpop and dec >= 0.5; keep.append(ok)
    print('gate %-8s hard members %6d (minpop %d), decisiveness %.2f -> %s' % (nm, nh, minpop, dec, 'MOTION' if ok else 'rejected'))
sites = [x for x, k in zip(sites, keep) if k]
if not sites:
    print('FACTORS %s 0 (no root-level site passes the gate: the product step adds nothing)' % tag); sys.exit(0)
n = len(sites); adj = np.eye(n, dtype=bool)
for i, j in itertools.combinations(range(n), 2):
    a, b = sites[i][1], sites[j][1]; ok = (a >= 0) & (b >= 0); m = nmi(a[ok] >= .5, b[ok] >= .5); r = np.corrcoef(a[ok], b[ok])[0, 1]
    dep = m > 0.05 or abs(r) > 0.2; adj[i, j] = adj[j, i] = dep
    print('  %-8s %-8s NMI %.3f r %+.2f %s' % (sites[i][0], sites[j][0], m, r, 'DEPENDENT' if dep else ''))
grp = -np.ones(n, int); g = 0
for i in range(n):
    if grp[i] >= 0: continue
    st = [i]; grp[i] = g
    while st:
        u = st.pop()
        for v in np.where(adj[u])[0]:
            if grp[v] < 0: grp[v] = g; st.append(v)
    g += 1
rep = [max(np.where(grp == k)[0], key=lambda i: margin.get(sites[i][0], 0)) for k in range(g)]
prod = 1 + sum(((sites[i][1] >= 0.5).astype(int) << k) for k, i in enumerate(rep))
cnt = np.bincount(prod)[1:]
print('motions %d: %s' % (g, '; '.join('{%s} -> %s' % (','.join(sites[i][0] for i in np.where(grp == k)[0]), sites[rep[k]][0]) for k in range(g))))
print('product: %d combinations, sizes %s' % (len(cnt), ' '.join(str(x) for x in cnt)))
for c, v in cols.items():
    ok = np.isfinite(v) & gen
    u = np.unique(v[ok])
    if 2 <= len(u) <= 20 and np.allclose(u, np.round(u)) and c not in ('junk',):
        print('  NMI with %s: product %.3f, leaves %.3f' % (c, nmi(prod[ok], v[ok]), nmi(lab[ok], v[ok])))
if 'labels_out' in args: np.savetxt(args['labels_out'], prod, fmt='%d')
print('FACTORS %s %d (root sites %d)' % (tag, g, n))
