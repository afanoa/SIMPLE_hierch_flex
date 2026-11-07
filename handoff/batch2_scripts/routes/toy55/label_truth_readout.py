#!/usr/bin/env python3
"""Truth readout of per-particle labels (one int per project row, 0 = none) on a synthetic benchmark: for 'arms' the
std ratio of each arm's angles within each label (relative to all particles); for 'tethered' the mean occupancy per site.
usage: label_truth_readout.py LABELS GT_PARTICLES TRUTH MODE"""
import sys, numpy as np
lab = np.loadtxt(sys.argv[1], dtype=int); gtf, trf, mode = sys.argv[2:5]
conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')}
tr = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
if mode.startswith('tethered'): tr = {c: v[0::4] for c, v in tr.items()}
d = len(next(iter(tr.values()))); N = max(len(lab), max(conf) + 1); T = np.full((N, d), np.nan)
for p, c in conf.items():
    if c in tr and p - 1 < N: T[p - 1] = tr[c]
lab = np.r_[lab, np.zeros(N - len(lab), int)]; allstd = np.nanstd(T, 0); ids = sorted(set(lab[lab > 0].tolist()))
print('%d labels; angle/occupancy std over all particles: %s' % (len(ids), ' '.join('%.2f' % v for v in allstd)))
best = []
for k in ids:
    t = T[lab == k]; t = t[~np.isnan(t[:, 0])]
    if len(t) < 30: continue
    if mode.startswith('arms'): r = np.nanstd(t, 0) / allstd; j = int(np.argmin(r)); best.append((r[j], k, len(t), j))
    else: m = t.mean(0); j = int(np.argmax(np.abs(m - 0.5))); best.append((1 - 2 * abs(m[j] - 0.5), k, len(t), j))
best.sort()
for v, k, n, j in best[:12]: print('  label %4d (%5d particles): best %s %d, score %.2f' % (k, n, 'arm' if mode.startswith('arms') else 'site', j // 2 + 1 if mode.startswith('arms') else j + 1, v))
print('summary: %d labels with >= 30 particles; median best score %.2f, best %.2f (arms: std ratio, 1 = no information; tethered: 1 - 2|occ - 0.5|, 0 = pure)' % (len(best), np.median([b[0] for b in best]) if best else 1, best[0][0] if best else 1))
