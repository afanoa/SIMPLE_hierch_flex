#!/usr/bin/env python3
"""orc_wrong.py PREP OUT [EVIDENCE]: stage B's wrong non-junk particles (their state's majority non-junk class is not theirs, as in
mest_post.py) against the right ones: pose error against the truth (in the stage-B frame), shift error, simulated quality (snrq) and
B-factor; the wrong share per pose-error and snrq bin. 'merged' = in their class's majority state, which another class owns (the
state count's limit); 'misassigned' = elsewhere. With EVIDENCE (stageB/prob_state_evidence.txt: iptcl label - ev_1..ev_K, the
last iteration's log-evidence) also the posteriors: confidently wrong = wrong with max posterior >= .99, and the mass they put
on a state their class owns."""
import os
import sys
import numpy as np

prep, fout = sys.argv[1:3]
fev = sys.argv[3] if len(sys.argv) > 3 else None
info = np.loadtxt(os.path.join(prep, 'ptcl_info.txt'), ndmin=2)
combo, junk, jtype = info[:, 0].astype(int), info[:, 1] > 0, info[:, 2].astype(int)
snrq, bfac, sB, ang, sh = info[:, 3], info[:, 4], info[:, 5].astype(int), info[:, 6], info[:, 7]
act = sB > 0
nj = act & ~junk
fo = open(fout, 'w')


def P(s=''):
    fo.write(s + '\n'); print(s, flush=True)


def auc(score, pos, neg):
    s = np.r_[score[pos], score[neg]]; r = np.argsort(np.argsort(s)) + 1.
    n1, n0 = pos.sum(), neg.sum()
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2.) / max(n1 * n0, 1))


maj = {}
for u in np.unique(sB[nj]):
    v = combo[nj & (sB == u)]
    maj[u] = np.bincount(v).argmax()
right = np.array([maj.get(u, -1) for u in sB]) == combo
own = {c: int(np.bincount(sB[nj & (combo == c)]).argmax()) for c in np.unique(combo[nj])}
ownstate = np.array([own.get(c, -1) for c in combo])
wrong = nj & ~right
merged = wrong & (sB == ownstate)
misas = wrong & ~merged
P('# %s: non-junk %d; stage-B wrong %.3f (merged %.3f, misassigned %.3f); state owners %s'
  % (prep, nj.sum(), wrong[nj].mean(), merged[nj].mean(), misas[nj].mean(), {int(u): int(c) for u, c in maj.items()}))
P('# classes whose majority state another class owns: %s' % sorted(int(c) for c in own if maj.get(own[c]) != c))


def row(name, m):
    if not m.any():
        return '%-22s n      0' % name
    return '%-22s n %6d | angle median %5.2f mean %5.2f p90 %6.2f >10 %.3f >30 %.3f | shift median %.2f p90 %.2f | snrq median %.3f mean %.3f | bfac median %.1f' % (
        name, m.sum(), np.median(ang[m]), ang[m].mean(), np.percentile(ang[m], 90), (ang[m] > 10).mean(), (ang[m] > 30).mean(),
        np.median(sh[m]), np.percentile(sh[m], 90), np.median(snrq[m]), snrq[m].mean(), np.median(bfac[m]))


P(row('right (non-junk)', nj & right))
P(row('wrong: merged', merged))
P(row('wrong: misassigned', misas))
P(row('junk', act & junk))
P('AUROC for misassigned vs right: angle error %.3f, shift error %.3f, low snrq %.3f, B-factor %.3f'
  % (auc(ang, misas, nj & right), auc(sh, misas, nj & right), auc(-snrq, misas, nj & right), auc(bfac, misas, nj & right)))
P()
P('# wrong share of non-junk by stage-B pose error (deg): n, merged, misassigned')
for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 30), (30, 181)):
    m = nj & (ang >= lo) & (ang < hi)
    if m.any():
        P('  [%3d,%3d) n %6d  merged %.3f  misassigned %.3f' % (lo, hi, m.sum(), merged[m].mean(), misas[m].mean()))
P('# wrong share of non-junk by snrq quintile')
q = np.percentile(snrq[nj], [0, 20, 40, 60, 80, 100])
for lo, hi in zip(q[:-1], q[1:]):
    m = nj & (snrq >= lo) & (snrq <= hi)
    P('  snrq [%.3f,%.3f] n %6d  merged %.3f  misassigned %.3f  angle median %.2f' % (lo, hi, m.sum(), merged[m].mean(), misas[m].mean(), np.median(ang[m])))
P('# per class: n, right, merged, misassigned; angle median right / misassigned')
for c in np.unique(combo[nj]):
    m = nj & (combo == c)
    P('  class %3d n %6d right %.3f merged %.3f misassigned %.3f | angle %.2f / %.2f' % (
        c, m.sum(), right[m].mean(), merged[m].mean(), misas[m].mean(), np.median(ang[m & right]) if (m & right).any() else np.nan,
        np.median(ang[m & misas]) if (m & misas).any() else np.nan))
if fev:
    d = np.loadtxt(fev, comments='#', ndmin=2)
    ipt = d[:, 0].astype(int); ev = np.full((len(sB), d.shape[1] - 3), np.nan); labe = np.zeros(len(sB), int)
    ok = (ipt >= 1) & (ipt <= len(sB)); ev[ipt[ok] - 1] = d[ok, 3:]; labe[ipt[ok] - 1] = d[ok, 1].astype(int)
    P()
    P('# posteriors from %s: evidence label = project label on %.4f of the active rows' % (fev, (labe[act] == sB[act]).mean()))
    e = np.where(np.isnan(ev) | (ev < -1e8), -np.inf, ev); mx = e.max(1, keepdims=True); mx[~np.isfinite(mx)] = 0.
    w = np.exp(e - mx); s = w.sum(1, keepdims=True); s[s == 0] = 1.; post = w / s
    mp = post.max(1)
    K = post.shape[1]
    okst = np.array([[maj.get(k + 1, -1) == combo[i] for k in range(K)] for i in range(len(sB))])
    mass = (post * okst).sum(1)
    conf = mp >= .99
    P('non-junk with max posterior >= .99: %.3f; of them wrong %.3f; wrong share among all non-junk %.3f'
      % (conf[nj].mean(), wrong[nj & conf].mean(), wrong[nj].mean()))
    for name, m in (('wrong', wrong), ('merged', merged), ('misassigned', misas)):
        if m.any():
            P('%-12s n %6d: max posterior median %.3f, >= .99 %.3f; mass on a state their class owns median %.3f mean %.3f'
              % (name, m.sum(), np.median(mp[m]), conf[m].mean(), np.median(mass[m]), mass[m].mean()))
    P(row('confidently right', nj & right & conf))
    P(row('confidently wrong', wrong & conf))
    P(row('  conf. misassigned', misas & conf))
    P(row('unconfident (<.99)', nj & ~conf))
    P('junk: max posterior mean %.3f, >= .99 %.3f' % (mp[act & junk].mean(), conf[act & junk].mean()))
fo.close()
