#!/usr/bin/env python3
"""Reference distance matrix (relative RMS of the 10 A conformer maps) among chosen conformers. usage: truth_dist_subset.py NPZ a,b,c,..."""
import sys, numpy as np
z = np.load(sys.argv[1]); ids = list(z['ids']); D = z['D']; S = [int(v) for v in sys.argv[2].split(':')]
print('      ' + ' '.join('%5d' % k for k in S))
for i in S: print('%5d ' % i + ' '.join('%5.2f' % D[ids.index(i), ids.index(j)] for j in S))
