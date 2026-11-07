#!/bin/bash
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Where does base80 lose the closed state? m33's two conformer-recovering closed leaves (r1111221_s1 / _s2) and all its closed-arm >= 0.6
# leaves, followed by project row into base80: their base80 leaves (count, purity, closed-arm) and every base80 node on their way down
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import os, numpy as np
P = '/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/'; S = '/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/'
def tree(run):
    T = P + run + '/1_hetero_tree'; L = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int); nm = [l.split()[0] for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
    hd, R = None, {}
    for l in open(S + run + '/fa_leaves/maps_table.txt'):
        w = l.split()
        if l.startswith('#'): hd = w[1:]; continue
        if hd and w: R[w[0]] = dict(zip(hd, w))
    return T, L, nm, R
mT, mL, mn, mR = tree('full692_m33_s1'); bT, bL, bn, bR = tree('full692_base80_s1'); NM = len(mL)
idx = lambda nm, s: nm.index(s) + 1
C = np.isin(mL, [idx(mn, 'r1111221_s1'), idx(mn, 'r1111221_s2')]); C2 = np.isin(mL, [idx(mn, k) for k, v in mR.items() if k in mn and float(v['closed_arm']) >= 0.6])
print('m33: the two closed leaves hold %d rows; all %d closed-arm >= 0.6 leaves hold %d' % (C.sum(), sum(1 for k, v in mR.items() if k in mn and float(v['closed_arm']) >= 0.6), C2.sum()))
for nmC, X in (('m33 closed leaves (conf18 / 19 matches)', C), ('m33 closed-arm >= 0.6 leaves', C2)):
    print('\n== %s -> base80 leaves (rows of them, share of the base80 leaf, its closed-arm / s1 / paired d):' % nmC)
    c = np.bincount(bL[X], minlength=len(bn) + 1); tot = np.bincount(bL, minlength=len(bn) + 1)
    for k in np.argsort(-c)[:8]:
        if c[k] == 0: break
        n = bn[k - 1] if k > 0 else '(no leaf)'; r = bR.get(n, {})
        print('   %-16s %6d of its %6d (%.2f)   %.0f %% of them   closed-arm %s s1 %s d %s' % (n, c[k], tot[k], c[k] / max(tot[k], 1), 100. * c[k] / X.sum(), r.get('closed_arm', '-'), r.get('s1_A', '-'), r.get('paired_d', '-')))
def nodes(T):
    out = {}
    for d in os.listdir(T):
        if d.startswith('node_r') and os.path.exists('%s/%s/rows.txt' % (T, d)): out[d[5:]] = np.loadtxt('%s/%s/rows.txt' % (T, d), dtype=int) - 1
    return out
for run, T, X in (('base80', bT, C), ('m33', mT, C)):
    ND = nodes(T); print('\n== %s nodes holding the m33 closed-leaf rows (%d): node, its rows, closed rows in it, their share of the node, their recall' % (run, X.sum()))
    cur = 'r'
    while True:
        r = ND[cur]; c = X[r].sum(); print('   %-10s %7d  %6d  %.3f  %.3f' % (cur, len(r), c, c / len(r), c / X.sum()))
        kids = [k for k in ND if k.startswith(cur) and len(k) == len(cur) + 1]
        if not kids: break
        for k in sorted(kids): rk = ND[k]; ck = X[rk].sum(); print('      child %-10s %7d  %6d  share %.3f  recall %.3f' % (k, len(rk), ck, ck / len(rk), ck / X.sum()))
        cur = max(kids, key=lambda k: X[ND[k]].sum())
PYEOF
