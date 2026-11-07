#!/usr/bin/env python3
"""Per-half meta for a gold-standard f692e run: every per-row array of META (first dimension = meta rows) restricted to
the half's rows in the half project's order (META_ROWS, 1-based meta rows); other entries copied.
usage: make_half_meta.py META META_ROWS OUT.npz"""
import sys, numpy as np
meta, mrf, out = sys.argv[1:4]; M = np.load(meta, allow_pickle=True); mr = np.loadtxt(mrf).astype(int) - 1; n = len(M['flags'])
D = {k: (M[k][mr] if (M[k].ndim >= 1 and M[k].shape[0] == n) else M[k]) for k in M.files}
np.savez(out, **D); print('%s: %d of %d rows; per-row keys %s; flags==1 %d' % (out, len(mr), n, [k for k in M.files if M[k].ndim >= 1 and M[k].shape[0] == n], int((D['flags'] == 1).sum())))
