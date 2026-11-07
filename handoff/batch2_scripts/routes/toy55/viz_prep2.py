#!/usr/bin/env python3
"""UMAP of a node's flex latent on its NPC leading view-centred principal components at their natural scale (PC scores, not
standardised per dimension), NSUB subsample; per particle the f692rf depth-2 branch (1 = r1, 2 = r21, 3 = r22, 0 = outside), the old
f692r closed flag (states 15+16) and the old f692r family (1 / 2). Writes OUT.npz (xy, branch, oldclosed, fam, evr = explained shares).
usage: viz_prep2.py FLEX OUT.npz [--npc 6 --nsub 80000 --oris ptcl3D_dump.txt]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('flex'); ap.add_argument('out'); ap.add_argument('--npc', type=int, default=6); ap.add_argument('--nsub', type=int, default=80000)
ap.add_argument('--oris', default='/mnt/beegfs/elmlund/afan/rnasep_klem/cheapkids_f692rf_r22/ptcl3D_parent.txt'); a = ap.parse_args()
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; R = K_ + '/rflex_f692rf'
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); cl = (old == 15) | (old == 16); NR = len(old)
fam = np.zeros(NR, int)
for f_, p in ((1, K_ + '/f692r/coordtree_fam1/units.txt'), (2, K_ + '/f692r/coordtree_fam2/units.txt')):
    try: u = np.loadtxt(p, dtype=int); fam[(u > 0) & (fam == 0)] = f_
    except Exception: pass
br = np.zeros(NR, int)
for b, nm in ((1, 'r1'), (2, 'r21'), (3, 'r22')): br[np.loadtxt(R + '/rows_%s.txt' % nm) > 0] = b
F = np.loadtxt(a.flex, comments='#', ndmin=2); pidx = F[:, 0].astype(int); Z = F[:, 5:]; N = len(pidx)
e1 = np.zeros(NR); e2 = np.zeros(NR); i = 0
for ln in open(a.oris):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); e1[i] = float(kv['e1']); e2[i] = float(kv['e2']); i += 1
rng = np.random.default_rng(0); E1, E2 = np.radians(e1[pidx - 1]), np.radians(e2[pidx - 1]); v = np.c_[np.sin(E2) * np.cos(E1), np.sin(E2) * np.sin(E1), np.cos(E2)]
nb = max(20, N // 1300); C = v[rng.choice(N, nb, replace=False)].copy()
for it in range(30):
    lb = np.argmax(v @ C.T, 1)
    for k in range(nb):
        mm = lb == k; C[k] = v[mm].mean(0) / np.linalg.norm(v[mm].mean(0)) if mm.any() else v[rng.integers(N)]
lb = np.argmax(v @ C.T, 1); _, bi = np.unique(lb, return_inverse=True); cnt = np.bincount(bi)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=bi.max() + 1) for j in range(Z.shape[1])], 1) / cnt[:, None]; Zc = Z - bm[bi]; Zc -= Zc.mean(0)
w, V = np.linalg.eigh(np.cov(Zc.T)); o = np.argsort(-w); w, V = w[o], V[:, o]; S = Zc @ V[:, :a.npc]
sub = np.sort(rng.choice(N, min(a.nsub, N), replace=False)); rows = pidx[sub] - 1
import umap
Y = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(S[sub])
np.savez_compressed(a.out, xy=Y.astype(np.float32), branch=br[rows], oldclosed=cl[rows], fam=fam[rows], evr=(w / w.sum())[:10], pc=S[sub, :3].astype(np.float32))
print('%s: %d particles, %d-dim latent, PC shares %s; UMAP of %d on %d PCs; branches %s; old closed %d'
      % (a.flex, N, Z.shape[1], ' '.join('%.2f' % x for x in (w / w.sum())[:8]), len(sub), a.npc, np.bincount(br[rows], minlength=4).tolist(), cl[rows].sum()))
