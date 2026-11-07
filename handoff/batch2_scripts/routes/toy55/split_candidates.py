#!/usr/bin/env python3
"""Candidate groups of a klem board for a data-derived split: the board's grid cells clustered into K groups (spectral
embedding of the cell-deviation correlation graph + k-means, as expand_families.py), particles by their argmax cell
(free-component particles 0). Writes a states file (one label per project row) and clusters.txt (state ncells mass n).
usage: split_candidates.py BOARD_DIR NROWS K OUT_DIR [MINSIZE]  (K lowered until every group has >= MINSIZE particles)"""
import sys, os, glob, numpy as np
bd, nrows, K, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; minsize = int(sys.argv[5]) if len(sys.argv) > 5 else 0; os.makedirs(out, exist_ok=True)
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
fs = sorted(glob.glob(bd + '/cells/cell_*.mrc')); V = np.array([read_mrc(f) for f in fs]); C = len(V); n = V.shape[1]; Dv = V - V.mean(0)
zz, yy, xx = np.indices((n, n, n)); band = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) >= 1
X = np.array([np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(d)))[band] for d in Dv]); X /= np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); S = np.real(X @ np.conj(X).T)
W = np.maximum(S, 0); np.fill_diagonal(W, 0); dg = W.sum(1) + 1e-9; w_, v_ = np.linalg.eigh(np.eye(C) - W / np.sqrt(np.outer(dg, dg)))
L = np.loadtxt(bd + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); on = am <= C
while True:
    emb = v_[:, 1:max(K, 2)]; best = None
    for seed in range(20):   # k-means on the spectral embedding, best of 20 starts
        rng = np.random.default_rng(seed); Cc = emb[rng.choice(C, K, replace=False)]
        for it in range(100): lab = ((emb[:, None, :] - Cc[None]) ** 2).sum(2).argmin(1); Cc = np.array([emb[lab == j].mean(0) if (lab == j).any() else Cc[j] for j in range(K)])
        J = ((emb - Cc[lab]) ** 2).sum()
        if best is None or J < best[0]: best = (J, lab.copy())
    lab = best[1]; st = np.zeros(nrows, int); st[pind[on] - 1] = lab[am[on] - 1] + 1
    if K <= 2 or min(int((st == j + 1).sum()) for j in range(K)) >= minsize: break
    K -= 1
np.savetxt(out + '/states.txt', st, fmt='%d')
with open(out + '/clusters.txt', 'w') as f:
    f.write('# state ncells mass nselected cells\n')
    for j in range(K): f.write('%d %d %d %d %s\n' % (j + 1, (lab == j).sum(), (st == j + 1).sum(), (st == j + 1).sum(), ' '.join(str(c + 1) for c in np.where(lab == j)[0])))
print('%d cells -> %d groups: cells %s, particles %s (free %d)' % (C, K, np.bincount(lab, minlength=K).tolist(), [int((st == j + 1).sum()) for j in range(K)], int((~on).sum())))
