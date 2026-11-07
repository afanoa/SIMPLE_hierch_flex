#!/usr/bin/env python3
"""usage: arm_score.py TREE_DIR TRUTH   (IgM-type toys: per-particle arm angles a<k>_th / a<k>_ph and a junk flag per project row)
Arm recovery of a hetero_tree run against the per-particle truth, over the genuine (truth non-junk) particles in leaves:
- eta^2 of every arm angle explained by leaf membership (between-leaf over total variance) and by unit membership (final_maps/
  units_index.txt; equal to the leaves' when every unit holds one leaf);
- span: the leaves' mean angles (leaves holding >= 1 % of the genuine particles in leaves) against the truth 10th / 90th percentiles of
  that angle over all genuine particles: the leaf-mean range over the P10-P90 range, and whether the lowest / highest leaf mean reaches
  P10 / P90;
- junk: the truth junk in junk-majority leaves, in the leaves the tree names junk and unassigned; whether every junk-majority leaf is
  named junk; the junk share inside the genuine-majority leaves;
- the delivered maps' half-map FSC line of final_maps/summary.txt.
ARMS lines (the scorecard reads them), then the per-angle table."""
import os, re, sys, numpy as np
T, TF = sys.argv[1].rstrip('/'), sys.argv[2]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; Tt = np.genfromtxt(TF, comments='#'); rr = Tt[:, 0].astype(int) - 1
L = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int); NM = len(L)
def col(c): v = np.full(NM, np.nan); v[rr] = Tt[:, hd.index(c)]; return v
junk = col('junk') > 0; ang = [h for h in hd if re.match(r'a\d+_(th|ph)$', h)]; gen = ~junk & ~np.isnan(col('junk'))
LV = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]; nm = {k + 1: w[0] for k, w in enumerate(LV)}; kd = {k + 1: w[1] for k, w in enumerate(LV)}
U = np.zeros(NM, int); ui = T + '/final_maps/units_index.txt'; idx = {w[0]: k + 1 for k, w in enumerate(LV)}
if os.path.exists(ui):
    for l in open(ui):
        w = l.split()
        if l.startswith('#') or len(w) < 2: continue
        for lf in re.split(r'[,+]', w[1]):
            if lf in idx: U[L == idx[lf]] = int(w[0])
def eta2(v, lab, sel):
    x, k = v[sel], lab[sel]; m = x.mean(); c = np.bincount(k); s = np.bincount(k, x); h = c > 0
    return (c[h] * (s[h] / c[h] - m) ** 2).sum() / ((x - m) ** 2).sum()
inl = gen & (L > 0); big = [k for k in nm if (inl & (L == k)).sum() >= 0.01 * inl.sum()]; rows = []
for a in ang:
    v = col(a); s = inl & np.isfinite(v); p10, p90 = np.nanpercentile(v[gen], [10, 90]); mu = [v[s & (L == k)].mean() for k in big]
    rows.append((a, eta2(v, L, s), eta2(v, U, s & (U > 0)) if (U > 0).any() else np.nan, p10, p90, min(mu), max(mu), (max(mu) - min(mu)) / (p90 - p10), min(mu) <= p10, max(mu) >= p90))
e, eu, cov = np.array([r[1] for r in rows]), np.array([r[2] for r in rows]), np.array([r[7] for r in rows]); nreach = sum(r[8] + r[9] for r in rows)
lv = [k for k in nm if (L == k).any()]; js = {k: junk[L == k].mean() for k in lv}; jm = [k for k in lv if js[k] > 0.5]; jn = [k for k in lv if kd[k] == 'junk']; gm = [k for k in lv if k not in jm]
J = junk.sum(); sz = sorted(((L == k).sum(), nm[k]) for k in lv)[-1]
print('ARMS leaves %d (%d named junk), units %d; %s particles in leaves; largest leaf %s (%s)' % (len(lv), len(jn), len(set(U[U > 0])), '{:,}'.format(int((L > 0).sum())), sz[1], '{:,}'.format(sz[0])))
print('ARMS arm angles over the genuine particles in leaves: eta^2 by leaf mean %.3f (min %.3f, max %.3f), by unit mean %.3f' % (e.mean(), e.min(), e.max(), np.nanmean(eu)))
print('ARMS span: leaf-mean range over the truth P10-P90 range, mean %.2f (min %.2f); %d of %d angle ends reached (leaves >= 1 %% of the genuine)' % (cov.mean(), cov.min(), nreach, 2 * len(rows)))
print('ARMS junk: %s truth junk; %.0f %% in %d junk-majority leaves (%s), %.0f %% in the %d leaves named junk, %.0f %% unassigned; junk-majority leaves named junk: %d of %d; junk share in the genuine-majority leaves %.3f'
      % ('{:,}'.format(int(J)), 100. * junk[np.isin(L, jm)].sum() / J, len(jm), ', '.join('%s %s %.2f%s' % (nm[k], '{:,}'.format(int((L == k).sum())), js[k], '' if kd[k] == 'junk' else ' unflagged') for k in jm),
         100. * junk[np.isin(L, jn)].sum() / J, len(jn), 100. * junk[L == 0].sum() / J, sum(kd[k] == 'junk' for k in jm), len(jm), junk[np.isin(L, gm)].sum() / max(np.isin(L, gm).sum(), 1)))
fl = [k for k in jn if k not in jm]
if fl: print('ARMS named junk but genuine-majority: %s' % ', '.join('%s %s %.2f' % (nm[k], '{:,}'.format(int((L == k).sum())), js[k]) for k in fl))
sm = T + '/final_maps/summary.txt'
if os.path.exists(sm): [print('ARMS maps: ' + l.strip()) for l in open(sm) if 'FSC' in l]
print('# angle eta2_leaf eta2_unit truth_P10 truth_P90 leaf_mean_min leaf_mean_max span reach_P10 reach_P90')
for r in rows: print('%-6s %.3f %.3f %+7.2f %+7.2f %+7.2f %+7.2f %.2f %d %d' % r)
