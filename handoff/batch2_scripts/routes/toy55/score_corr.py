#!/usr/bin/env python3
"""usage: score_corr.py NAME=FILE ..   (scores: one value per row, < 0 outside the level)
The Spearman matrix of the scores over the rows where all of them are >= 0, to show which candidate motions are the same one."""
import sys, numpy as np
from scipy.stats import spearmanr
nm = [a.split('=')[0] for a in sys.argv[1:]]; X = np.column_stack([np.loadtxt(a.split('=', 1)[1]) for a in sys.argv[1:]]); ok = (X >= 0).all(1)
R = spearmanr(X[ok])[0] if X.shape[1] > 2 else np.array([[1., spearmanr(X[ok, 0], X[ok, 1])[0]], [spearmanr(X[ok, 0], X[ok, 1])[0], 1.]])
print('Spearman over %d rows:' % ok.sum()); print('%-10s ' % '' + ' '.join('%8s' % x for x in nm))
for i, a in enumerate(nm): print('%-10s ' % a + ' '.join('%+8.3f' % R[i, j] for j in range(len(nm))))
