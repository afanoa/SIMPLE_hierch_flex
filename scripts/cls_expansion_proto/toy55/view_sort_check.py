#!/usr/bin/env python3
"""Are the leaf states view-sorted? Per leaf: the projection-direction histogram (100 spherical k-means bins of the whole
set) of its particles vs the whole set: max bin ratio, KL divergence, and the fraction of the leaf in its 10 most
over-represented bins. usage: view_sort_check.py RUN META_NPZ [--maxleaves 40]"""
import sys, os, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('meta'); ap.add_argument('--maxleaves', type=int, default=40); a = ap.parse_args()
M = np.load(a.meta); e1, e2 = np.radians(M['e1']), np.radians(M['e2']); v = np.c_[np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)]; N = len(v)
rng = np.random.default_rng(0); C = v[rng.choice(N, 100, replace=False)].copy()
for it in range(30):
    lab = np.argmax(v @ C.T, 1)
    for k in range(100):
        m = lab == k
        if m.any(): C[k] = v[m].mean(0); C[k] /= np.linalg.norm(C[k])
lab = np.argmax(v @ C.T, 1); H0 = np.bincount(lab, minlength=100) / N
rows = []
for cf in sorted(glob.glob(a.run + '/fam?/region*/level1/states_clusters.txt'))[:a.maxleaves]:
    st = np.loadtxt(cf, dtype=int); name = re.search(r'(fam\d/region[^/]+)/', cf).group(1)
    for s in sorted(set(st[st > 0])):
        m = np.where(st == s)[0]; H = np.bincount(lab[m], minlength=100) / len(m); r = H / np.maximum(H0, 1e-9); kl = float(np.sum(H[H > 0] * np.log(H[H > 0] / H0[H > 0]))); top = np.sort(r)[::-1]
        rows.append((name, s, len(m), kl, float(r.max()), float(H[np.argsort(-r)[:10]].sum())))
print('leaf  particles  KL(view | whole)  max bin ratio  share in its 10 most over-represented bins (uniform = 0.10)')
for r in rows: print('%-26s s%02d  %6d  %.3f  %5.2f  %.2f' % r)
kls = np.array([r[3] for r in rows]); print('summary: %d leaves; KL median %.3f (a random subset of the same size: ~%.3f); max bin ratio median %.2f' % (len(rows), np.median(kls), 100 / (2 * np.median([r[2] for r in rows])), np.median([r[4] for r in rows])))
