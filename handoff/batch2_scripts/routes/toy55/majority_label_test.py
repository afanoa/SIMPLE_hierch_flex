#!/usr/bin/env python3
"""usage: majority_label_test.py TRUTH TCOL TABLE SITE: a two-pole site's label sizes from the port's size/site_size.txt EPMs: each pole's own
top-EPM vs sizing the MAJORITY pole as the complement of the minority pole's EPM (top (n_level - EPM_minority) by the majority pole's
posterior); precision / recall per pole of the pole's truth class (TCOL > 0 for the pole that tracks it, == 0 for the other)."""
import sys, os, numpy as np
TJ, TC, TAB, SITE = sys.argv[1:5]; hd = [l.split()[1:] for l in open(TJ) if l.startswith('#')][0]; T = np.loadtxt(TJ); c = lambda n: T[:, hd.index(n)]
junk = c('junk') > 0; y1 = (c(TC) > 0) & ~junk; y0 = (c(TC) == 0) & ~junk; NM = len(T)
w = [l.split() for l in open(TAB) if l.startswith(SITE + ' ')][0]; Q = np.loadtxt(w[3], ndmin=2); r = np.loadtxt(w[4], dtype=int).ravel()
rows = np.zeros(NM, bool); rows[r - 1] = True if not (len(r) == NM and r.max() <= 1) else False
if len(r) == NM and r.max() <= 1: rows = r > 0
S = os.path.join(os.path.dirname(os.path.dirname(w[2])), '..', 'node_r', 'fsite_' + SITE, 'size', 'site_size.txt')
print('site', SITE, 'level rows', rows.sum(), '| site_size.txt candidates:', [p for p in (S,) if os.path.exists(p)])
def topk(P, k): o = np.argsort(-np.where(rows & (P >= 0), P, -np.inf)); s = np.zeros(NM, bool); s[o[:k]] = True; return s
def pr(s, y): return (s & y).sum() / max((s & ~junk).sum(), 1), (s & y).sum() / max((rows & y).sum(), 1)
for j in range(Q.shape[1]):
    P = Q[:, j]; ok = rows & (P >= 0)
    for y, nm in ((y1, TC + '=1'), (y0, TC + '=0')):
        o = np.argsort(-np.where(ok, P, -np.inf)); top = y[o[:2500]].mean()
        print('  pole %d posterior: top-2.5k purity for %s %.3f' % (j + 1, nm, top))
EPM = [int(v) for v in sys.argv[5].split(',')] if len(sys.argv) > 5 else None
if EPM:
    n = rows.sum(); mino = int(np.argmin(EPM)); majo = 1 - mino; ymin, ymaj = (y1, y0) if Q[:, mino][y1 & rows].mean() > Q[:, mino][y0 & rows].mean() else (y0, y1)
    TL = np.loadtxt(sys.argv[6], dtype=int) if len(sys.argv) > 6 else np.ones(NM, int); scr = rows & (TL > 0)
    pmin, pmaj = Q[:, mino], Q[:, majo]; smin = topk(pmin, EPM[mino]) & scr; arg = scr & (pmaj > pmin)
    for lab, s_, y in (('two-way argmax, junk-screened: pole %d (majority)' % (majo + 1), arg, ymaj), ('two-way argmax, junk-screened: pole %d (minority)' % (mino + 1), scr & ~arg, ymin),
                       ('min top-EPM + argmax of the rest: majority', scr & ~smin & (pmaj > pmin), ymaj),
                       ('complement within junk-screened: majority (%d)' % (scr.sum() - (smin).sum()), scr & ~smin, ymaj)):
        p_, r_ = pr(s_, y); print('  %-48s precision %.3f recall %.3f (%d labelled)' % (lab, p_, r_, s_.sum()))
    for lab, s, y in (('minority pole %d, own top-EPM %d' % (mino + 1, EPM[mino]), topk(Q[:, mino], EPM[mino]), ymin),
                      ('majority pole %d, own top-EPM %d' % (majo + 1, EPM[majo]), topk(Q[:, majo], EPM[majo]), ymaj),
                      ('majority pole %d, complement n - EPM_min = %d' % (majo + 1, n - EPM[mino]), topk(Q[:, majo], n - EPM[mino]), ymaj)):
        p_, r_ = pr(s, y); print('  %-48s precision %.3f recall %.3f' % (lab, p_, r_))
