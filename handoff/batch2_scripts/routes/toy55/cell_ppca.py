#!/usr/bin/env python3
"""Cell-weighted PPCA of a finished klem board = the node's locally-linear chart from ALL its particles, no new
reconstruction: the particle's expected deviation is its (row-normalised) responsibility-weighted combination of the
cell deviations, so the particle covariance in deviation space is D^H C D with C = Rn^T Rn (ngrid x ngrid); PCs from
the SVD of C^(1/2) (D - mu_w); particle coordinates z = Rn (D - mu_w) V. Cross-half reproducibility per PC from D0/D1.
Controls on the same particles: the unweighted cell PCA (the current level-2 axis) and a flex_pca embedding (--flexdir).
Truth readout (--gt, --states, --axis): Spearman(PC, truth axis), 15-NN conformer accuracy (balanced), closed-vs-rest AUC.
usage: cell_ppca.py BOARD OUT [--npc 4] [--gt gt_particles.txt --states 1,10,.. --axis a1,.. --closed 18,51] [--flexdir DIR]"""
import os, sys, argparse, numpy as np
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('out'); ap.add_argument('--npc', type=int, default=4); ap.add_argument('--gt', default=''); ap.add_argument('--states', default=''); ap.add_argument('--axis', default='')
ap.add_argument('--closed', default='18,51'); ap.add_argument('--flexdir', default=''); ap.add_argument('--minmass', type=float, default=0.5); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
H = np.load(a.board + '/klem_halves.npz'); D0, D1 = H['D0'].astype(complex), H['D1'].astype(complex); ngrid = len(D0); R = np.load(a.board + '/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int)
gm = r[:, :ngrid].sum(1); ok = gm >= a.minmass; Rn = r[ok, :ngrid] / gm[ok, None]; pinds = pinds[ok]; n = len(pinds); D = 0.5 * (D0 + D1)
def chart(Dm, weighted):
    w = Rn.sum(0); mu = (w @ Dm) / w.sum() if weighted else Dm.mean(0); Dc = Dm - mu
    if weighted:
        C = Rn.T @ Rn; ev, U = np.linalg.eigh(C); A = (U * np.sqrt(np.maximum(ev, 0))) @ U.T; B = A @ Dc
    else: B = Dc
    _, S, Vh = np.linalg.svd(B, full_matrices=False); V = Vh[:a.npc].conj().T; T = Dc @ V
    rep = [np.real(np.vdot((D0 - mu) @ V[:, j], (D1 - mu) @ V[:, j])) / max(np.linalg.norm((D0 - mu) @ V[:, j]) * np.linalg.norm((D1 - mu) @ V[:, j]), 1e-30) for j in range(a.npc)]
    return np.real(Rn @ T), S[:a.npc] / S[0], rep
zw, sw, repw = chart(D, True); zu, su, repu = chart(D, False)
np.savetxt(a.out + '/cell_ppca_z.txt', np.c_[pinds, zw], fmt='%d' + ' %.5f' * a.npc, header='pind z1..z%d (cell-weighted PPCA)' % a.npc)
print('board %s: %d particles on the grid (%d cells); weighted PPCA singular values %s, cross-half reproducibility of the PCs %s' % (a.board.rstrip('/').split('/')[-1], n, ngrid, np.array2string(sw, precision=2), np.array2string(np.array(repw), precision=2)))
print('  unweighted cell PCA (current level-2 axis): singular values %s, reproducibility %s' % (np.array2string(su, precision=2), np.array2string(np.array(repu), precision=2)))
if not a.gt: sys.exit(0)
gt = {}
for l in open(a.gt):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = (int(w[1]), int(w[2]))
states = [int(v) for v in a.states.split(',')]; axis = dict(zip(states, [float(v) for v in a.axis.split(',')])); closed = set(int(v) for v in a.closed.split(','))
conf = np.array([gt.get(p, (0, 0))[0] for p in pinds]); m = conf > 0; ax = np.array([axis.get(c, 0.) for c in conf]); cl = np.array([c in closed for c in conf])
def readout(z, tag):
    zz = (z - z.mean(0)) / np.maximum(z.std(0), 1e-12); sp = [abs(spearmanr(zz[m, j], ax[m])[0]) for j in range(zz.shape[1])]
    tree = cKDTree(zz[m]); _, nb = tree.query(zz[m], k=16); pred = np.array([np.bincount(conf[m][row[1:]]).argmax() for row in nb]); acc = (pred == conf[m]).mean()
    bal = np.mean([(pred[conf[m] == c] == c).mean() for c in states if (conf[m] == c).any()]); chance = max(np.bincount(conf[m])) / m.sum()
    auc = max(abs(0.5 - (np.argsort(np.argsort(zz[m, j]))[cl[m]].mean() / m.sum())) * 2 for j in range(zz.shape[1]))   # rank AUC distance from chance, best PC
    print('  %-28s Spearman(PC, truth axis) per PC %s | 15-NN accuracy %.3f (balanced %.3f, chance %.3f) | closed-vs-rest |AUC-0.5|x2 %.2f' % (tag, np.array2string(np.array(sp), precision=2), acc, bal, chance, auc))
readout(zw, 'cell-weighted PPCA'); readout(zu, 'unweighted cell PCA')
if a.flexdir:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from latent_windows import read_embedding
    _, fp, fz, _ = read_embedding(a.flexdir + '/flex_pca_embedding.bin'); pos = {int(p): i for i, p in enumerate(fp)}; idx = np.array([pos.get(int(p), -1) for p in pinds]); keep = idx >= 0
    zf = np.zeros((n, min(a.npc, fz.shape[1]))); zf[keep] = fz[idx[keep], :zf.shape[1]]; readout(zf, 'flex_pca latent (global)')
