#!/usr/bin/env python3
"""Pick a realistic known-structure cut from the 55-conformer landscape: the most populated conformer of each reference
clade (average linkage of the true-map distances cut at --cut), a within-clade partner at a moderate distance for the two
largest clades, and one near-identical partner (should merge). usage: pick_confs.py TRUTH_NPZ GT_A GT_B [--cut 0.40]"""
import sys, argparse, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
ap = argparse.ArgumentParser(); ap.add_argument('npz'); ap.add_argument('gta'); ap.add_argument('gtb'); ap.add_argument('--cut', type=float, default=0.40)
ap.add_argument('--mid', default='0.33,0.45'); ap.add_argument('--near', type=float, default=0.20); a = ap.parse_args()
z = np.load(a.npz); D = z['D']; ids = list(z['ids']); cl = fcluster(linkage(squareform(D, checks=False), 'average'), t=a.cut, criterion='distance')
cnt = {h: np.bincount(np.loadtxt(f, comments='#')[:, 1].astype(int), minlength=max(ids) + 1) for h, f in (('A', a.gta), ('B', a.gtb))}
n = {k: min(cnt['A'][k], cnt['B'][k]) for k in ids}; ix = {k: i for i, k in enumerate(ids)}; lo, hi = [float(v) for v in a.mid.split(',')]
clades = sorted(set(cl), key=lambda c: -sum(n[k] for k in ids if cl[ix[k]] == c)); pick = []
for ci, c in enumerate(clades):
    mem = sorted([k for k in ids if cl[ix[k]] == c], key=lambda k: -n[k]); rep = mem[0]; pick.append(rep)
    print('clade %d: %d conformers, %d particles per half (min of halves); representative c%d (%d per half)' % (c, len(mem), sum(n[k] for k in mem), rep, n[rep]))
    if ci < 2:
        mids = [k for k in mem[1:] if lo <= D[ix[rep], ix[k]] <= hi]
        if mids: p = max(mids, key=lambda k: n[k]); pick.append(p); print('   moderate partner c%d at distance %.2f (%d per half)' % (p, D[ix[rep], ix[p]], n[p]))
    if ci == 0:
        nears = [k for k in mem[1:] if D[ix[rep], ix[k]] < a.near and k not in pick]
        if nears: p = max(nears, key=lambda k: n[k]); pick.append(p); print('   near-identical partner c%d at distance %.2f (%d per half)' % (p, D[ix[rep], ix[p]], n[p]))
print('picked:', ':'.join(str(k) for k in pick), '| particles per half:', sum(n[k] for k in pick))
print('pairwise distances:'); print('      ' + ' '.join('%6s' % ('c%d' % k) for k in pick))
for k in pick: print('%6s' % ('c%d' % k) + ' '.join('%6.2f' % D[ix[k], ix[j]] for j in pick))
