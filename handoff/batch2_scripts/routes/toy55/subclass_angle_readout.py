#!/usr/bin/env python3
"""Per-angle readout of 2D sub-class labels on a hinged-arm synthetic set: for every label with >= MINP particles the std
of each arm's in-plane (theta) and out-of-plane (phi) angle relative to all conformer particles. Per arm: the best theta
and best phi ratio over labels; over labels: the median of each label's best theta and best phi ratio and how many labels
pin some arm's theta / phi below 0.5 and 0.3 (1 = no information). usage: subclass_angle_readout.py LABELS GT TRUTH [MINP 30]"""
import sys, numpy as np
labf, gtf, trf = sys.argv[1:4]; minp = int(sys.argv[4]) if len(sys.argv) > 4 else 30
lab = np.loadtxt(labf, dtype=int, ndmin=1); G = np.loadtxt(gtf, comments='#'); conf = np.zeros(int(G[:, 0].max()) + 1, int); conf[G[:, 0].astype(int)] = G[:, 1].astype(int)
ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
narm = len(next(iter(ang.values()))) // 2; c = conf[1:len(lab) + 1]; ok = c > 0
A = np.full((len(lab), 2 * narm), np.nan); A[ok] = np.array([ang[k] for k in c[ok]]); sall = np.nanstd(A[ok], 0)
R, sizes = [], []
for l in sorted(set(lab[lab > 0].tolist())):
    rows = np.where((lab == l) & ok)[0]
    if len(rows) >= minp: R.append(np.nanstd(A[rows], 0) / sall); sizes.append(len(rows))
R = np.array(R); th, ph = R[:, 0::2], R[:, 1::2]; bt, bp = th.min(1), ph.min(1)
print('%d labels with >= %d conformer particles (median size %d); %d arms' % (len(R), minp, int(np.median(sizes)), narm))
for k in range(narm):
    i, j = int(np.argmin(th[:, k])), int(np.argmin(ph[:, k]))
    print('  arm %2d: best in-plane theta ratio %.2f (label of %d), best out-of-plane phi ratio %.2f (label of %d)' % (k + 1, th[i, k], sizes[i], ph[j, k], sizes[j]))
print('per label, best arm: theta median %.2f, phi median %.2f; labels pinning a theta < 0.5: %d, < 0.3: %d; a phi < 0.5: %d, < 0.3: %d'
      % (np.median(bt), np.median(bp), (bt < 0.5).sum(), (bt < 0.3).sum(), (bp < 0.5).sum(), (bp < 0.3).sum()))
