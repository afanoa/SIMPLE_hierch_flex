#!/usr/bin/env python3
"""usage: qtoy_score.py TREE_DIR TRUTH   (qtoy: run_qtoy_sim.sbatch's truth_qtoy.txt; TREE_DIR = .../1_hetero_tree)
Truth-scored readouts of a hetero_tree run on qtoy (states plus a per-particle quality spread independent of them):
- per state (the arm-path conformers, the closed conformer, the holo, blob and broken junk): the leaves where it is the plurality, its
  recall (share of its particles in them) and their precision (its share of their particles);
- arm: eta^2 of the truth s1 by leaf over the arm-path particles; the leaf means of s1 (leaves with >= 1 % of them) against the path's
  lowest / highest conformer;
- quality: eta^2 of the truth B and log q by leaf over the particles in leaves (about 0 when the leaves ignore quality);
- per split (the run's tree.log 'split group' lines, node split/split_children.txt): eta^2 between the children of B and log q (quality)
  and of s1 (arm-path particles), closed, holo and junk (state); the call is the largest, none under 0.02;
- junk: in junk-majority leaves, in leaves named junk, unassigned; the junk share of the other leaves.
QTOY lines (the scorecard reads them), then the per-state and per-split tables."""
import os, re, sys, numpy as np
T, TF = sys.argv[1].rstrip('/'), sys.argv[2]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; Tt = np.genfromtxt(TF, comments='#')
L = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int); NM = len(L); rw = Tt[:, hd.index('row')].astype(int) - 1; ok = rw < NM
def col(c): v = np.full(NM, np.nan); v[rw[ok]] = Tt[ok, hd.index(c)]; return v
conf, jt, s1, Bf, q = (col(c) for c in ('conformer', 'jtype', 's1', 'bfac', 'snrq')); junk, holo, closed = col('junk') > 0, col('holo') > 0, col('C_closed') > 0
st = np.where(junk, -np.nan_to_num(jt), np.nan_to_num(conf)).astype(int); arm = ~junk & ~holo & ~closed & np.isfinite(s1); lq = np.log(q)
LV = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]; nm = {k + 1: w[0] for k, w in enumerate(LV)}; kd = {k + 1: w[1] for k, w in enumerate(LV)}
lv = [k for k in nm if (L == k).any()]
def eta2(v, lab, sel):
    sel = sel & np.isfinite(v)
    if sel.sum() < 20: return np.nan
    x, k = v[sel], lab[sel]; m = x.mean(); c = np.bincount(k); s = np.bincount(k, x); h = c > 0; t = ((x - m) ** 2).sum()
    return (c[h] * (s[h] / c[h] - m) ** 2).sum() / t if t > 0 else np.nan
def sname(s):
    if s == -1: return 'blob'
    if s == -2: return 'broken'
    m = st == s; return 'conf%02d %s' % (s, 'closed' if closed[m].any() else 'holo' if holo[m].any() else 's1 %+.1f' % np.nanmean(s1[m]))
pl = {k: int(np.bincount(st[L == k] + 2).argmax()) - 2 for k in lv}; states = sorted(set(st.tolist()), key=lambda s: (s < 0, np.nanmean(s1[st == s]) if (st == s).any() and np.isfinite(s1[st == s]).any() else 99, s))
SR = []
for s in states:
    lk = [k for k in lv if pl[k] == s]; inl = np.isin(L, lk); n = (st == s).sum()
    SR.append((sname(s), n, len(lk), (inl & (st == s)).sum() / n, (inl & (st == s)).sum() / max(inl.sum(), 1), (L[st == s] == 0).mean()))
ia = arm & (L > 0); big = [k for k in lv if (ia & (L == k)).sum() >= 0.01 * ia.sum()]; mu = sorted((s1[ia & (L == k)].mean(), nm[k]) for k in big)
pth = (np.nanmin(s1[arm]), np.nanmax(s1[arm])); il = L > 0
jm = [k for k in lv if junk[L == k].mean() > 0.5]; jn = [k for k in lv if kd[k] == 'junk']; J = junk.sum(); gm = [k for k in lv if k not in jm]
mB = sorted((Bf[L == k].mean(), np.exp(lq[L == k].mean()), nm[k], (L == k).sum()) for k in lv)
cs = [x for x in SR if 'closed' in x[0]]; cb = max(lv, key=lambda k: (closed & (L == k)).sum())
g = ((cs[0][0], cs[0][3], cs[0][4], cs[0][2]) if cs else ('closed', 0., 0., 0)) + (nm[cb], '{:,}'.format(int((L == cb).sum())), closed[L == cb].mean(), 100. * (closed & (L == cb)).sum() / max(closed.sum(), 1))
print('QTOY CLOSED GATE: %s: recall %.2f / precision %.2f in %d leaves; most of it in %s (%s particles, closed share %.2f, %.0f %% of the closed)' % g)
print('QTOY leaves %d (%d named junk), %s particles in leaves, %s unassigned' % (len(lv), len(jn), '{:,}'.format(int(il.sum())), '{:,}'.format(int((~il).sum()))))
print('QTOY states (recall / precision through the leaves each one dominates, in N leaves): ' + ' | '.join('%s %s: %.2f / %.2f in %d' % (a, '{:,}'.format(int(n)), r, p, nl) for a, n, nl, r, p, u in SR))
print('QTOY arm: eta^2 of s1 by leaf %.3f over the arm-path particles; leaf means %+.1f (%s) .. %+.1f (%s) A vs the path %+.1f .. %+.1f A (%.0f %% of it)' % (eta2(s1, L, ia),
      mu[0][0], mu[0][1], mu[-1][0], mu[-1][1], pth[0], pth[1], 100. * (mu[-1][0] - mu[0][0]) / (pth[1] - pth[0])))
print('QTOY quality by leaf: eta^2 of B %.3f, of log q %.3f; leaf-mean B %.0f .. %.0f A^2, q %.2f .. %.2f (truth: B mean %.0f, q median %.2f)' % (eta2(Bf, L, il), eta2(lq, L, il),
      mB[0][0], mB[-1][0], min(x[1] for x in mB), max(x[1] for x in mB), np.nanmean(Bf), np.nanmedian(q)))
print('QTOY junk: %s truth junk; %.0f %% in %d junk-majority leaves (%s), %.0f %% in the %d leaves named junk, %.0f %% unassigned; junk-majority leaves named junk %d of %d; junk share of the other leaves %.3f'
      % ('{:,}'.format(int(J)), 100. * junk[np.isin(L, jm)].sum() / J, len(jm), ', '.join('%s %s %.2f%s' % (nm[k], '{:,}'.format(int((L == k).sum())), junk[L == k].mean(), '' if kd[k] == 'junk' else ' unnamed') for k in jm),
         100. * junk[np.isin(L, jn)].sum() / J, len(jn), 100. * junk[~il].sum() / J, sum(kd[k] == 'junk' for k in jm), len(jm), junk[np.isin(L, gm)].sum() / max(np.isin(L, gm).sum(), 1)))
ga, gi = ~junk & ~il, ~junk & il
print('QTOY set aside (in no leaf): %s rows, %s truth junk and %s truth-genuine; the genuine set aside: B mean %.0f vs %.0f in leaves, q median %.2f vs %.2f%s'
      % ('{:,}'.format(int((~il).sum())), '{:,}'.format(int((junk & ~il).sum())), '{:,}'.format(int(ga.sum())), np.nanmean(Bf[ga]) if ga.any() else np.nan, np.nanmean(Bf[gi]),
         np.nanmedian(q[ga]) if ga.any() else np.nan, np.nanmedian(q[gi]), ('; by state: ' + ', '.join('%s %d' % (sname(s), (ga & (st == s)).sum()) for s in states if (ga & (st == s)).any())) if ga.any() else ''))
SP, TL = {}, os.path.join(os.path.dirname(T), 'tree.log')
for l in (open(TL, errors='replace') if os.path.exists(TL) else []):
    m = re.search(r'HETERO_TREE node (\S+): split group (\d+) \((\d+) particles; [^)]*\) -> (.*)$', l)
    if m: SP.setdefault(m.group(1), {})[int(m.group(2))] = m.group(4).strip().replace(' ', '_')
SL = []
for X in sorted(SP, key=lambda s: (len(s), s)):
    N = '%s/node_%s' % (T, X); f = N + '/split/split_children.txt'
    if not os.path.exists(f): continue
    ch = np.loadtxt(f, dtype=int)
    if len(ch) != NM: c2 = np.zeros(NM, int); c2[np.loadtxt(N + '/rows.txt', dtype=int) - 1] = ch; ch = c2
    sel = ch > 0; e = {'B': eta2(Bf, ch, sel), 'logq': eta2(lq, ch, sel), 's1': eta2(s1, ch, sel & arm), 'closed': eta2(closed.astype(float), ch, sel),
                       'holo': eta2(holo.astype(float), ch, sel), 'junk': eta2(junk.astype(float), ch, sel)}
    k = max(e, key=lambda k: e[k] if np.isfinite(e[k]) else -1); call = 'none' if not e[k] >= 0.02 else 'quality' if k in ('B', 'logq') else 'arm' if k == 's1' else k
    SL.append((X, int(sel.sum()), len(SP[X]), e, call))
cc = {}; [cc.__setitem__(x[4], cc.get(x[4], 0) + 1) for x in SL]
print('QTOY splits %d by truth (eta^2 between the children; the largest of quality B / log q, arm s1, closed, holo, junk): %s' % (len(SL), ', '.join('%s %d' % kv for kv in sorted(cc.items(), key=lambda kv: -kv[1]))))
qs = [x for x in SL if x[4] == 'quality']
if qs: print('QTOY quality-following splits: ' + ', '.join('%s (%s; B %.2f, log q %.2f vs best state %.2f)' % (x[0], '{:,}'.format(x[1]), x[3]['B'], x[3]['logq'], max(v for kk, v in x[3].items() if kk not in ('B', 'logq') and np.isfinite(v))) for x in qs))
print('# state particles leaves recall precision unassigned')
for a, n, nl, r, p, u in SR: print('%-22s %6d %3d %.3f %.3f %.3f' % (a, n, nl, r, p, u))
print('# split node particles groups eta2_B eta2_logq eta2_s1 eta2_closed eta2_holo eta2_junk call')
for X, n, g, e, call in SL: print('%-10s %6d %d %s %s' % (X, n, g, ' '.join('%.3f' % e[k] for k in ('B', 'logq', 's1', 'closed', 'holo', 'junk')), call))
print('# leaf particles mean_B mean_q')
for b, qq, n_, c in mB: print('%-16s %6d %6.1f %.3f' % (n_, c, b, qq))
