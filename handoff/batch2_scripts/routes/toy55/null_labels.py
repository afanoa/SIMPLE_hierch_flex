#!/usr/bin/env python3
"""No-heterogeneity null for the final pass: the particles of ONE conformer (the most abundant open one unless given) split
into K random groups per half, every other row 0. usage: null_labels.py GT_A GT_B OUT_A OUT_B K [CONF]"""
import sys, numpy as np
ga, gb, oa, ob, K = sys.argv[1:6]; K = int(K); ca, cb = [np.loadtxt(g, comments='#') for g in (ga, gb)]
conf = int(sys.argv[6]) if len(sys.argv) > 6 else None
if conf is None:
    cnt = np.bincount(ca[(ca[:, 1] > 0) & (ca[:, 2] == 0), 1].astype(int)); conf = int(cnt.argmax())
for g, o, seed in ((ca, oa, 1), (cb, ob, 2)):
    lab = np.zeros(len(g), int); r = np.where(g[:, 1].astype(int) == conf)[0]; lab[r] = np.random.default_rng(seed).integers(1, K + 1, len(r))
    np.savetxt(o, lab, fmt='%d'); print('conformer %d: %d rows -> %d random groups %s' % (conf, len(r), K, np.bincount(lab[r])[1:].tolist()))
