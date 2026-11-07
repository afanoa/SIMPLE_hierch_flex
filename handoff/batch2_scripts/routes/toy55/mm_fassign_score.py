#!/usr/bin/env python3
"""Score a re-assignment (mm_fassign.sbatch states.txt, one state per project row, 0 = excluded) on the multi-motion toy, genuine rows:
 factor mode (table=comb_table.txt): each state's per-component bits; every truth factor is matched to its best component (truth scores
   only); per factor the accuracy of the bit, the two pooled groups' map contrast; combination NMI of the states;
 cell mode (no table): states are cells (leaves); per factor the cell-majority accuracy and the best cell pair's contrast; combination NMI.
 Map contrast = |frac(state 2 in A) - frac(state 2 in B)| x sqrt(nA nB / (nA + nB)) between the two delivered maps.
usage: mm_fassign_score.py STATES truth=.. factors=.. [table=comb_table.txt] [tag=..]"""
import sys, numpy as np
args = dict(a.split('=', 1) if '=' in a else ('states', a) for a in sys.argv[1:]); tag = args.get('tag', args['states'])
t = np.loadtxt(args['truth'], ndmin=2); rows = t[:, 0].astype(int); combo = t[:, 1].astype(int); junk = t[:, 2] > 0
F = np.loadtxt(args['factors'], ndmin=2); nf = F.shape[1] - 3; fmap = {int(x[0]): x[1:1 + nf].astype(int) for x in F}
fs = np.zeros((len(rows), nf), int)
for i in np.where(~junk)[0]: fs[i] = fmap[combo[i]]
st = np.loadtxt(args['states'], dtype=int, ndmin=1)[rows - 1]; g = ~junk & (st > 0)
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; p = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(p, (a, b), 1); p /= p.sum()
    pa, pb = p.sum(1), p.sum(0); nz = p > 0; mi = (p[nz] * np.log(p[nz] / (pa[:, None] * pb[None])[nz])).sum()
    return mi / max(np.sqrt(-(pa[pa > 0] * np.log(pa[pa > 0])).sum() * -(pb[pb > 0] * np.log(pb[pb > 0])).sum()), 1e-12)
def contrast(A, B, s2): return abs(s2[A].mean() - s2[B].mean()) * np.sqrt(A.sum() * B.sum() / (A.sum() + B.sum())) if A.sum() and B.sum() else 0.
out = []
if 'table' in args:
    tab = np.loadtxt(args['table'], dtype=int, ndmin=2); bits = {r_[0]: r_[1:] for r_ in tab}; r = tab.shape[1] - 1
    Bt = np.array([bits.get(s, np.zeros(r, int)) for s in st])
    for q in range(nf):
        s2 = fs[:, q] == 2; best = max((max(((Bt[g, j] == 2) == s2[g]).mean(), ((Bt[g, j] == 2) != s2[g]).mean()), j) for j in range(r)); j = best[1]
        A = g & (Bt[:, j] == 2); Bm = g & (Bt[:, j] == 1); out.append('f%d acc %.3f (c%d) contrast %.0f (n %d/%d)' % (q + 1, best[0], j + 1, contrast(A, Bm, s2), A.sum(), Bm.sum()))
else:
    cells = [c for c in np.unique(st[g])]
    for q in range(nf):
        s2 = fs[:, q] == 2; acc = sum(max(s2[g & (st == c)].sum(), (~s2[g & (st == c)]).sum()) for c in cells) / g.sum()
        best = max((contrast(g & (st == a), g & (st == b), s2), a, b) for i, a in enumerate(cells) for b in cells[i + 1:])
        out.append('f%d acc %.3f contrast %.0f (cells %d/%d)' % (q + 1, acc, *best))
print('%s: %d genuine rows in %d states; combination NMI %.3f; %s' % (tag, g.sum(), len(np.unique(st[g])), nmi(st[g], combo[g]), '; '.join(out)))
