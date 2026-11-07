#!/usr/bin/env python3
"""Control for mm_fassign.sbatch: the same one-pass re-assignment with the tree's OWN leaves as references (comb_NN.mrc = the leaf's
even + odd unfiltered half-map average, comb_states.txt = hetero_tree_labels.txt, 0 kept), so a factor re-assignment's gain over it
comes from the factor model, not from one more refinement pass. usage: mm_leaf_refs.py TREE_DIR OUT_DIR"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
T, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
L = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
def loc(p): return p if os.path.exists(p) else os.path.join(T, p[p.index('1_hetero_tree') + len('1_hetero_tree') + 1:])
for i, l in enumerate(L):
    v = 0.5 * (read_mrc(loc(l[5])).astype(np.float32) + read_mrc(loc(l[6])).astype(np.float32)); B = v.shape[0]
    write_mrc('%s/comb_%02d.mrc' % (out, i + 1), v, 160 * 1.464 / B)
lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1); np.savetxt(out + '/comb_states.txt', lab, fmt='%d')
print('%d leaf references, states %s' % (len(L), ' '.join('%d:%d' % (k, (lab == k).sum()) for k in range(0, len(L) + 1))))
