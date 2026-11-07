#!/usr/bin/env python3
"""usage: motion_bins5.py SCORE OUT_STATES EXCLUDE TREE_DIR [TRUTH COLUMN ..]   (EXCLUDE tree | truth | none)
motion_bins4.py's at-node bins without a fixed truth kind: 5 quantile bins of SCORE (one value per row, < 0 outside the level), ties
at random (seed 7), with the excluded rows at state 0 and the edges set on the remaining level rows:
- tree: rows the tree calls junk (no leaf, or a leaf TREE_DIR/hetero_tree_leaves.txt names 'junk');
- truth: the truth junk (an oracle, needs TRUTH);
- none: nothing is excluded (TREE_DIR is not read; pass -).
With TRUTH (a sites truth table with a junk column) it prints per bin the truth junk share and, among genuine rows, the share of each
COLUMN > 0 ('COLUMN 0.123', the form phaseE_maptests.py reads as a truth key)."""
import sys, numpy as np
SC, OUT, EX, D = sys.argv[1:5]; TF = sys.argv[5] if len(sys.argv) > 5 else None; COLS = sys.argv[6:]
sc = np.loadtxt(SC); NM = len(sc); lev = sc >= 0; junk = np.zeros(NM, bool); comp = {}
if TF:
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; M = np.genfromtxt(TF, comments='#'); tr = M[:, 0].astype(int) - 1
    assert tr.max() + 1 == NM, 'truth rows %d != score rows %d' % (tr.max() + 1, NM)
    def col(c): v = np.zeros(NM); v[tr] = M[:, hd.index(c)]; return v
    junk = col('junk') > 0; comp = {c: (col(c) > 0).astype(float) for c in COLS}
if EX == 'tree':
    L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int); assert len(L) == NM, 'tree rows %d != score rows %d' % (len(L), NM)
    LV = [l.split() for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]; out = (L == 0) | np.isin(L, [k + 1 for k, w in enumerate(LV) if w[1] == 'junk'])
else: out = junk.copy() if EX == 'truth' else np.zeros(NM, bool)
ref = lev & ~out; s = sc + np.random.default_rng(7).uniform(0., 1e-9, NM); sg = np.sort(s[ref]); n = len(sg)
cuts = [sg[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]; st = np.where(ref, np.searchsorted(np.array(cuts), s, side='left') + 1, 0); np.savetxt(OUT, st, fmt='%d')
print('%s: %d rows, %d in the level, %d excluded (%s), %d binned' % (SC, NM, lev.sum(), (lev & out).sum(), EX, n))
for k in range(1, 6):
    m = st == k; g = m & ~junk
    txt = (', junk share %.3f' % ((m & junk).sum() / max(m.sum(), 1)) if TF else '') + ''.join(', %s %.3f' % (c, v[g].mean()) for c, v in comp.items())
    print('bin %d: rows %6d, score %.4f-%.4f%s' % (k, m.sum(), sc[m].min(), sc[m].max(), txt))
