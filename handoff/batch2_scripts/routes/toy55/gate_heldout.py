#!/usr/bin/env python3
"""Held-out check of the site-product gate (>= minpop hard members and decisiveness = hard / soft pole size >= 0.5) on trees not used to
pick it. Per root-level kept site: the site layer's own numbers (factorial/sites.txt), the gate's call, and a truth class fixed in advance:
 qtoy2   REAL if pole NMI with conformer (genuine rows) >= .05 or |rho(posterior, closed_arm)| >= .3, and pole junk share < .5
 igm     REAL if max |rho(posterior, arm angle)| >= .3 and pole junk share < .5
 detlim  REAL if pole F1 with the rare conformer >= .3 or pole NMI with conformer >= .05
 mic3    REAL if pole NMI with the protein state >= .05
 null / qtoynull / micnull   every site FAKE (no structure)
 ole     no truth (reported only)
A disagreement is a real site the gate rejects (MISS) or a fake site it passes (FALSE PASS). Sites whose numbers equal a site seen while
choosing .5 are marked 'seen'. usage: gate_heldout.py RUNS_FILE(name path kind) ROOT [seen=FILE(name n z hard per line)]"""
import sys, os, re, numpy as np
from scipy.stats import spearmanr
runs, root = sys.argv[1], sys.argv[2]; args = dict(a.split('=', 1) for a in sys.argv[3:])
TR = {'qtoy2': root + '/../ctl/truth/qtoy2.txt', 'qtoynull': root + '/../ctl/truth/qtoynull.txt', 'igm': root + '/../igm/arms_truth.txt',
      'micnull': root + '/../ctl/truth/micnull.txt', 'mic3': root + '/truth/mic3.txt'}
seen = set()
if 'seen' in args:
    for l in open(args['seen']): w = l.split(); seen.add((w[1], round(float(w[2])), round(float(w[3]), 1), round(float(w[4]), 3)))
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; p = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(p, (a, b), 1); p /= p.sum()
    x, y = p.sum(1), p.sum(0); nz = p > 0; mi = (p[nz] * np.log(p[nz] / np.outer(x, y)[nz])).sum(); return mi / max(np.sqrt(-(x[x > 0] * np.log(x[x > 0])).sum() * -(y[y > 0] * np.log(y[y > 0])).sum()), 1e-12)
def truth(kind, nm):
    f = TR.get(kind) or (root + '/truth/%s.txt' % nm if kind == 'detlim' else None)
    if not f or not os.path.exists(f): return None
    hdr = [l for l in open(f) if l.startswith('#')][0].lstrip('#').split(); t = np.loadtxt(f, comments='#', ndmin=2); return {h: t[:, i] for i, h in enumerate(hdr)}
print('%-16s %-8s %7s %7s %6s %6s %-9s %-6s %s' % ('run', 'site', 'soft', 'hard', 'dec', 'z', 'gate', 'truth', 'evidence'))
tally = {}
for l in open(runs):
    nm, p, kind = l.split(); R = root + '/' + p; T = R + '/1_hetero_tree'
    if not os.path.exists(T + '/factorial/table.txt'): continue
    rootsites = [x.split() for x in open(T + '/factorial/table.txt') if not x.startswith('#') and x.split()[1] == 'r']
    if not rootsites: continue
    mf = ''.join(open(f).read() for f in (R + '/run_manifest.txt', R + '/manifest.txt') if os.path.exists(f)); mp = re.search(r'minpop=(\d+)', mf); minpop = int(mp.group(1)) if mp else 900
    lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1); N = len(lab)
    fr = np.loadtxt(R + '/fixture_rows.txt', dtype=int, ndmin=1) if os.path.exists(R + '/fixture_rows.txt') and len(lab) == sum(1 for _ in open(R + '/fixture_rows.txt')) else None
    tr = truth(kind, nm); st = {}
    for x in open(T + '/factorial/sites.txt'):
        if x.startswith('#') or '|' not in x: continue
        w = x.split('|'); fz = lambda v: 1e4 if '*' in v else float(v); st[w[0].split()[0]] = (fz(w[2].split()[0]), fz(w[3].split()[0]), fz(w[4].split()[0]))   # ***** = overflowed print
    for s in rootsites:
        name = s[0]; soft, z, hs = st.get(name, (0, 0, 0)); nh = hs * N; dec = nh / soft if soft else 0.
        gate = 'PASS' if nh >= minpop and dec >= .5 else 'reject'
        pf = s[3]; pf = pf if os.path.exists(pf) else T + '/factorial/post/' + os.path.basename(pf); post = np.loadtxt(pf, ndmin=2)[:, 0]
        cls, ev = 'n/a', ''
        if kind in ('null', 'qtoynull', 'micnull'): cls, ev = 'FAKE', 'null fixture'
        elif tr is not None:
            # project row of every labelled row -> truth row
            prow = fr if fr is not None else np.arange(1, N + 1); pos = {int(r): i for i, r in enumerate(tr['row'])}
            ti = np.array([pos.get(int(r), -1) for r in prow]); ok = (ti >= 0) & (post >= 0); ti = ti[ok]; pp = post[ok]; pole = pp >= .5
            junk = tr['junk'][ti] > 0 if 'junk' in tr else np.zeros(len(ti), bool); js = junk[pole].mean() if pole.any() else 0.; g = ~junk
            if kind == 'qtoy2':
                a = nmi(pole[g], tr['conformer'][ti][g]) if pole[g].any() else 0.; ca = tr['closed_arm'][ti]; fa = g & np.isfinite(ca); r = abs(spearmanr(pp[fa], ca[fa])[0])
                cls = 'REAL' if (a >= .05 or r >= .3) and js < .5 else 'FAKE'; ev = 'conf NMI %.3f, |rho arm| %.2f, pole junk %.2f, snrq rho %+.2f' % (a, r, js, spearmanr(pp, tr['snrq'][ti])[0])
            elif kind == 'igm':
                rs = {c: abs(spearmanr(pp[g], tr[c][ti][g])[0]) for c in tr if c.startswith('a')}; c = max(rs, key=rs.get)
                cls = 'REAL' if rs[c] >= .3 and js < .5 else 'FAKE'; ev = 'best %s |rho| %.2f, pole junk %.2f' % (c, rs[c], js)
            elif kind == 'detlim':
                rare = tr['conformer'][ti] == 18; f1 = 2. * (pole & rare).sum() / max(pole.sum() + rare.sum(), 1); a = nmi(pole[g], tr['conformer'][ti][g]) if pole[g].any() else 0.
                cls = 'REAL' if f1 >= .3 or a >= .05 else 'FAKE'; ev = 'rare F1 %.2f (rare %d), conf NMI %.3f, pole junk %.2f, snrq rho %+.2f' % (f1, rare.sum(), a, js, spearmanr(pp, tr['snrq'][ti])[0])
            elif kind == 'mic3':
                a = nmi(pole, tr['state'][ti]) if pole.any() else 0.; cls = 'REAL' if a >= .05 else 'FAKE'
                ev = 'state NMI %.3f, micpc1 rho %+.2f, bfac rho %+.2f, snrq rho %+.2f' % (a, spearmanr(pp, tr['micpc1'][ti])[0], spearmanr(pp, tr['bfac'][ti])[0], spearmanr(pp, tr['snrq'][ti])[0])
        tag = 'seen' if (name, round(soft), round(z, 1), round(hs, 3)) in seen else ''
        flag = ''
        if cls == 'REAL' and gate == 'reject': flag = 'MISS'
        if cls == 'FAKE' and gate == 'PASS': flag = 'FALSE PASS'
        tally.setdefault((cls, gate, bool(tag)), 0); tally[(cls, gate, bool(tag))] += 1
        print('%-16s %-8s %7d %7d %6.2f %6.1f %-9s %-6s %s %s %s' % (nm, name, soft, nh, dec, z, gate, cls, ev, tag, flag))
print('TALLY (truth, gate, seen-before): ' + '; '.join('%s/%s%s %d' % (k[0], k[1], ' seen' if k[2] else '', v) for k, v in sorted(tally.items())))
