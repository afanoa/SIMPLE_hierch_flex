#!/usr/bin/env python3
"""Per-particle merged-group labels of a split test: candidate states.txt + the merge npz (names, labels) -> labels file.
usage: group_labels.py STATES MERGE_NPZ OUT"""
import sys, numpy as np
st = np.loadtxt(sys.argv[1]).astype(int); z = np.load(sys.argv[2]); g = np.zeros(len(st), int)
for nm, lb in zip(z['names'], z['labels']): g[st == int(str(nm).split('_s')[-1])] = lb
np.savetxt(sys.argv[3], g, fmt='%d'); print(sys.argv[3], np.bincount(g).tolist())
