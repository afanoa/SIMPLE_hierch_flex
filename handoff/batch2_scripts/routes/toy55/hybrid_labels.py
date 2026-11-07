#!/usr/bin/env python3
"""Tree + site layer as one partition: every tree leaf split by the site-layer classes of its rows (the common refinement of the two
labellings); a piece below MINN rows stays with the largest piece of its leaf. Writes OUT (per row class id, 0 = none) and prints the
pieces. usage: hybrid_labels.py TREE_LABELS SITE_LABELS[,SITE_LABELS2 ...] OUT [--minn 500]  (site labels: 1 / 2 per row, 0 = outside)"""
import sys, numpy as np
tl, sl, out = sys.argv[1:4]; minn = int(sys.argv[sys.argv.index('--minn') + 1]) if '--minn' in sys.argv else 500
L = np.loadtxt(tl, dtype=int); S = [np.loadtxt(f, dtype=int) for f in sl.split(',')]; key = L.astype(np.int64)
for s in S: key = key * 3 + np.clip(s, 0, 2)
key[L <= 0] = -1; out_ = np.zeros(len(L), int); nid = 0; npieces = 0
for leaf in np.unique(L[L > 0]):
    m = L == leaf; ks, cnt = np.unique(key[m], return_counts=True); big = ks[np.argmax(cnt)]; ids = {}
    for k, c in zip(ks, cnt):
        if c >= minn or k == big: nid += 1; ids[k] = nid; npieces += 1
    for k in ks: out_[m & (key == k)] = ids.get(k, ids[big])
np.savetxt(out, out_, fmt='%d'); print('%d tree leaves x %d site labellings -> %d classes (pieces >= %d rows)' % (len(np.unique(L[L > 0])), len(S), nid, minn))
