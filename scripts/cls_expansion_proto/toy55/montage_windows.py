#!/usr/bin/env python3
"""Montage: window maps of a latent-windows arm beside the truth conformer that dominates each window, projections at
60-10 A (three axes), plus the window-minus-mean difference against the truth-minus-consensus difference.
usage: montage_windows.py MAPDIR INFO CONS --truth f1,.. --states .. --out fig.png [--nw 8] [--smpd 1.464]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, Shells
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('mapdir'); ap.add_argument('info'); ap.add_argument('cons'); ap.add_argument('--truth', required=True); ap.add_argument('--states', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--nw', type=int, default=8); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--title', default='')
a = ap.parse_args(); states = [int(s) for s in a.states.split(',')]; T = dict(zip(states, a.truth.split(',')))
cons = read_mrc(a.cons); B = cons.shape[0]; S = Shells(B, a.smpd); msk = S.soft_mask(100. / a.smpd)
def bl(v): return S.lowpass(S.lowpass(v * msk, 10) - S.lowpass(v * msk, 60), 10) * 0 + (S.lowpass(v * msk, 10) - S.lowpass(v * msk, 60))
info = np.loadtxt(a.info, comments='#', ndmin=2); n = info[:, 1]; dom = info[:, 2].astype(int); share = info[:, 3]
M = [read_mrc('%s/recvol_state%02d.mrc' % (a.mapdir, w)) for w in range(1, a.nw + 1)]; Mm = sum(ni * m for ni, m in zip(n, M)) / n.sum()
fig, ax = plt.subplots(4, a.nw, figsize=(2.2 * a.nw, 9.2))
for j in range(a.nw):
    W = bl(M[j]); Tm = bl(read_mrc(T[dom[j]])) if dom[j] in T else np.zeros_like(W); D = bl(M[j] - Mm); Dt = bl(read_mrc(T[dom[j]]) - cons) if dom[j] in T else np.zeros_like(W)
    for i, (V, cmap, name) in enumerate(((W, 'gray_r', 'window'), (Tm, 'gray_r', 'truth %d' % dom[j]), (D, 'RdBu_r', 'window - mean'), (Dt, 'RdBu_r', 'truth - consensus'))):
        P = V.sum(axis=0); vmax = np.abs(P).max() if cmap == 'RdBu_r' else None
        ax[i, j].imshow(P, cmap=cmap, vmin=-vmax if vmax else None, vmax=vmax); ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
        if i == 0: ax[i, j].set_title('w%d n=%d\ndominant %d (%.2f)' % (j + 1, n[j], dom[j], share[j]), fontsize=8)
        if j == 0: ax[i, j].set_ylabel(name, fontsize=9)
fig.suptitle(a.title or '%s: window maps vs dominant truth conformer, 60-10 A projections along z' % a.mapdir, fontsize=10); plt.tight_layout(); plt.savefig(a.out, dpi=100); print('wrote', a.out)
