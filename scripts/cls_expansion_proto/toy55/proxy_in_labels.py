#!/usr/bin/env python3
"""Share of the earlier pure closed sub-state particles (proxy) in every label of a half-row labelling.
usage: proxy_in_labels.py PROXY_LABELS_META 1:3 META_ROWS LABELS [NAME]"""
import sys, numpy as np
pl, want, mrf, labf = sys.argv[1:5]; name = sys.argv[5] if len(sys.argv) > 5 else labf
P = np.isin(np.loadtxt(pl).astype(int), [int(v) for v in want.split(':')]); mr = np.loadtxt(mrf).astype(int) - 1; px = P[mr]; L = np.loadtxt(labf).astype(int)
print('%s: proxy %d in the selection of %d; per label: %s' % (name, px[L > 0].sum(), (L > 0).sum(), ' '.join('%d:%d/%d (%.1f%%)' % (k, px[L == k].sum(), (L == k).sum(), 100. * px[L == k].mean()) for k in range(1, L.max() + 1))))
