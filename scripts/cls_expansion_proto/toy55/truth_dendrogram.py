#!/usr/bin/env python3
"""Reference hierarchy of a synthetic set's conformers: every conformer map (docked, box 160) low-passed to LP A and masked,
pairwise RMS difference relative to the mean map's RMS, average-linkage clustering; prints the merge heights at which the
conformers join (top-down splits with populations) and saves the distance matrix. usage: truth_dendrogram.py MAPDIR GT OUT_NPZ [LP]"""
import sys, glob, os, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
mapdir, gtf, out = sys.argv[1:4]; lp = float(sys.argv[4]) if len(sys.argv) > 4 else 10.; smpd = 1.464
conf = np.loadtxt(gtf, comments='#')[:, 1].astype(int); pop = np.bincount(conf); ids = [k for k in range(1, len(pop)) if pop[k] > 0]
fs = {k: (sorted(glob.glob('%s/conf%02d_*_m160.mrc' % (mapdir, k))) or sorted(glob.glob('%s/conf%02d*.mrc' % (mapdir, k))))[0] for k in ids}
V0 = read_mrc(fs[ids[0]]); n = V0.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (lp / smpd / 1.18) ** 2)
x = np.arange(n) - n // 2; rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = rr < 0.45 * n
M = np.array([np.real(np.fft.ifftn(np.fft.fftn(read_mrc(fs[k]).astype(np.float64)) * G))[msk] for k in ids]); M = (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)
D = np.sqrt(((M[:, None, :] - M[None, :, :]) ** 2).mean(2)); np.savez(out, ids=np.array(ids), D=D, pop=pop[ids])
Z = linkage(squareform(D, checks=False), 'average'); root = to_tree(Z)
def show(nd, depth=0, maxdepth=6):
    leaves = nd.pre_order(lambda x: ids[x.id]); p = sum(pop[k] for k in leaves)
    lab = ('%d conformers, %d particles' % (len(leaves), p)) if len(leaves) > 4 else ' '.join('%d(%d)' % (k, pop[k]) for k in leaves)
    print('%s%s  [join height %.3f]' % ('  ' * depth, lab, nd.dist))
    if not nd.is_leaf() and depth < maxdepth and len(leaves) > 1: show(nd.get_left(), depth + 1); show(nd.get_right(), depth + 1)
print('%d conformers; pairwise relative RMS difference at %.0f A: median %.3f, 10th pct %.3f, max %.3f' % (len(ids), lp, np.median(D[np.triu_indices(len(ids), 1)]), np.percentile(D[np.triu_indices(len(ids), 1)], 10), D.max()))
show(root)
