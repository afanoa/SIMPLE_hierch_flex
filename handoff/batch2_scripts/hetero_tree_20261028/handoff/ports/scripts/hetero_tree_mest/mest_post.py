#!/usr/bin/env python3
"""mest_post.py OUT.txt: stage-B posterior sharpness, evidence gaps, calibration and the soft-versus-hard membership purity of the
gate trees' nodes (mu_* trees: their stageB/prob_state_evidence.txt is refine3D's per-particle, per-state log-evidence of the last
iteration), the numbers behind 'the weights are one-hot'."""
import sys, os
import numpy as np

G = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates'
TQ = '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt'
CREF = '/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt'
BINS = [0, 0.9, 0.99, 0.999, 1.0000001]


def read_evidence(f, n):
    d = np.loadtxt(f, comments='#', ndmin=2); ipt = d[:, 0].astype(int); lab = np.zeros(n, int); ev = np.full((n, d.shape[1] - 3), np.nan)
    ok = (ipt >= 1) & (ipt <= n); lab[ipt[ok] - 1] = d[ok, 1].astype(int); ev[ipt[ok] - 1] = d[ok, 3:]
    return lab, ev


def posteriors(ev):
    e = np.where(np.isnan(ev) | (ev < -1e8), -np.inf, ev); mx = e.max(1, keepdims=True); mx[~np.isfinite(mx)] = 0.
    w = np.exp(e - mx); s = w.sum(1, keepdims=True); s[s == 0] = 1.
    return w / s


out = open(sys.argv[1], 'w')
P_ = lambda s: (out.write(s + '\n'), print(s, flush=True))
hdr = open(TQ).readline().lstrip('#').split(); tq = np.loadtxt(TQ, ndmin=2)
combo_all = tq[:, hdr.index('combo')].astype(int); junk_all = tq[:, hdr.index('junk')].astype(int)
cref = np.loadtxt(CREF, dtype=int)
fxr = np.loadtxt(os.path.join(G, 'mu_r131', 'fixture_rows.txt'), dtype=int)
P_('# stage-B posteriors (softmax of the per-state log-evidence) per node: n, K, mean max posterior, shares >= .9 / .99 / .999, gap = best - second evidence (nats) p10/p50/p90')
P_('# qtoy2: calibration = accuracy of the argmax state (its majority combo) per max-posterior bin [<.9 | .9-.99 | .99-.999 | >.999], then the misassigned share and the posterior mass they put on a state of their own combo; purity hard vs soft-weighted')
P_('# r131: reference closed rows per state (enrichment over the node share) and the closed rows\' mean max posterior')
for tree, nodes in (('mu_qtoy2', ['r', 'r1', 'r2', 'r3', 'r4', 'r5']), ('mu_r131', ['r', 'r1', 'r2']), ('mu_qtoynull', ['r']), ('mu_olenull', ['r'])):
    for nd in nodes:
        d = os.path.join(G, tree, '1_hetero_tree', 'node_' + nd); rows = np.loadtxt(os.path.join(d, 'rows.txt'), dtype=int, ndmin=1); n = len(rows)
        lab, ev = read_evidence(os.path.join(d, 'stageB', 'prob_state_evidence.txt'), n); P = posteriors(ev); m = lab > 0
        K = P.shape[1]; mp = P[m].max(1); srt = -np.sort(-np.where(np.isnan(ev[m]), -1e9, ev[m]), 1); gap = srt[:, 0] - srt[:, 1]
        line = '%-14s %6d %2d meanmax %.3f  >.9 %.3f  >.99 %.3f  >.999 %.3f | gap %6.1f %6.1f %6.1f' % (tree[3:] + '/' + nd, m.sum(), K, mp.mean(), (mp >= .9).mean(), (mp >= .99).mean(), (mp >= .999).mean(), *np.percentile(gap, [10, 50, 90]))
        if tree == 'mu_qtoy2':
            combo = combo_all[rows - 1]; junk = junk_all[rows - 1]; L = lab[m]; C = combo[m]; J = junk[m] > 0; Pm = P[m]
            maj = {k: np.bincount(C[(L == k) & ~J]).argmax() for k in np.unique(L)}
            right = np.array([maj[k] for k in L]) == C
            acc = []
            for lo, hi in zip(BINS[:-1], BINS[1:]):
                b = (mp >= lo) & (mp < hi) & ~J; acc.append('%.3f[%d]' % (right[b].mean() if b.sum() else float('nan'), b.sum()))
            wrong = ~right & ~J
            ok_state = np.array([[maj.get(k + 1, -1) == C[i] for k in range(K)] for i in range(len(C))])
            hp = (Pm * ok_state).sum(1)[wrong]
            hard_pur = sum(((L == k) & (C == maj[k]) & ~J).sum() for k in maj) / (~J).sum()
            soft_pur = sum((Pm[~J, k - 1] * (C[~J] == maj[k])).sum() for k in maj) / Pm[~J].sum()
            line += '\n      calibration %s | wrong %.3f of non-junk, their mass on a right state: median %.3f mean %.3f | junk share %.3f, junk meanmax %.3f | purity hard %.3f soft %.3f' % (
                ' '.join(acc), wrong.sum() / (~J).sum(), np.median(hp) if len(hp) else 0, hp.mean() if len(hp) else 0, J.mean(), mp[J].mean(), hard_pur, soft_pur)
        if tree == 'mu_r131':
            full = fxr[rows - 1]; isref = np.isin(full, cref); L = lab[m]; R = isref[m]; sh = R.mean()
            enr = ['s%d: n %d ref %d (x%.2f)' % (k, (L == k).sum(), (R & (L == k)).sum(), ((R & (L == k)).sum() / max((L == k).sum(), 1)) / sh) for k in np.unique(L)]
            line += '\n      closed refs %d (%.3f): %s; ref meanmax %.3f others %.3f' % (R.sum(), sh, ' '.join(enr), mp[R].mean(), mp[~R].mean())
        P_(line)
out.close()
