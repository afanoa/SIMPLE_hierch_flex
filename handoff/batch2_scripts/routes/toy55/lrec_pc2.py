#!/usr/bin/env python3
# lrec_pc.py with leaves as k-means clusters in the whitened space of the top NPC principal components of the centred score
# vectors (a family with contamination or two motions: the motion is not PC1 alone); NPC 1 = quantile bins of PC1 as lrec_pc.py.
# usage: lrec_pc2.py OUTDIR --tabs 'dir/dist_part_neigh_*.dat' --nstates K --nrows N --k KL --npc 2 [--prev coord.txt] [--gt gt.txt]
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--tabs', required=True); ap.add_argument('--nstates', type=int, required=True)
ap.add_argument('--nrows', type=int, required=True); ap.add_argument('--k', type=int, default=8); ap.add_argument('--npc', type=int, default=2)
ap.add_argument('--prev'); ap.add_argument('--gt'); a = ap.parse_args()
CAND = np.dtype([('iref', '<i4'), ('inpl', '<i4'), ('dist', '<f4'), ('x', '<f4'), ('y', '<f4'), ('has_sh', '<i4')])
P, Ds = [], []
for f in sorted(glob.glob(a.tabs)):
    sz = os.path.getsize(f); raw = open(f, 'rb').read(); nrefs, nptcls, nnz, nchunks = np.frombuffer(raw, '<i8', 4)
    nspace = nrefs // a.nstates; assert nspace * a.nstates == nrefs, (nrefs, a.nstates)
    pinds = np.frombuffer(raw, '<i4', nptcls, 32); pos = sz - (28 * nnz + 4 * nchunks); D = np.full((nptcls, a.nstates), np.inf, np.float32)
    for c in range(nchunks):
        n = int(np.frombuffer(raw, '<i4', 1, pos)[0]); pos += 4; idx = np.frombuffer(raw, '<i4', n, pos) - 1; pos += 4 * n
        cand = np.frombuffer(raw, CAND, n, pos); pos += 24 * n; st = (cand['iref'] - 1) // nspace
        ok = (cand['inpl'] > 0) & (st >= 0) & (st < a.nstates); np.minimum.at(D, (idx[ok], st[ok]), cand['dist'][ok])
    assert pos == sz, (f, pos, sz); P.append(pinds.copy()); Ds.append(D)
pind = np.concatenate(P); D = np.concatenate(Ds).astype(float); full = np.isfinite(D).all(1); pind, D = pind[full], D[full]
Dc = D - D.mean(1, keepdims=True); Dc -= Dc.mean(0); U, S, Vt = np.linalg.svd(Dc, full_matrices=False); Z = U[:, :a.npc] * S[:a.npc]
print('tables: %d particles; PCA variance fractions: %s' % (len(pind), ' '.join('%.3f' % v for v in S[:5] ** 2 / (S ** 2).sum())))
if a.npc == 1:
    q = np.quantile(Z[:, 0], np.linspace(0, 1, a.k + 1)[1:-1]); b = np.searchsorted(q, Z[:, 0]) + 1
else:
    X = Z / Z.std(0); rng = np.random.default_rng(0); best = None
    for rep in range(8):
        C = X[rng.choice(len(X), a.k, replace=False)]
        for it in range(100):
            lab = ((X[:, None] - C[None]) ** 2).sum(-1).argmin(1); C = np.array([X[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(a.k)])
        sse = ((X - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
    b = best[1] + 1
os.makedirs(a.out, exist_ok=True); lab = np.zeros(a.nrows, int); lab[pind - 1] = b; np.savetxt(os.path.join(a.out, 'bins.txt'), lab, fmt='%d')
for j in range(a.npc):
    pc = np.full(a.nrows, np.nan); pc[pind - 1] = Z[:, j]; np.savetxt(os.path.join(a.out, 'pc%d.txt' % (j + 1)), pc, fmt='%.6g')
np.savetxt(os.path.join(a.out, 'coord.txt'), np.loadtxt(os.path.join(a.out, 'pc1.txt')), fmt='%.6g')
if a.gt:
    gt = np.loadtxt(a.gt, comments='#', usecols=(1,), dtype=int)[pind - 1]; sw = (gt >= 200) & (gt <= 230); x = gt - 200.
    W = sum(((x[sw & (b == v)] - x[sw & (b == v)].mean()) ** 2).sum() for v in np.unique(b[sw])); e = 1 - W / ((x[sw] - x[sw].mean()) ** 2).sum()
    w = np.linalg.lstsq(np.c_[Z[sw], np.ones(sw.sum())], x[sw], rcond=None)[0][:a.npc]
    rows = sorted([(x[sw & (b == v)].mean(), (b == v).sum(), sw[b == v].mean()) for v in np.unique(b) if (sw & (b == v)).any()])
    print('TRUTH: swing share %.3f; corr(PC1..%d, swing) %s, best direction in them %.3f; %d leaves swing eta2 %.3f; leaves (mean swing/n/swing share): %s'
          % (sw.mean(), a.npc, ' '.join('%.3f' % abs(np.corrcoef(Z[sw, j], x[sw])[0, 1]) for j in range(a.npc)), np.corrcoef(Z[sw] @ w, x[sw])[0, 1], a.k, e,
             ' '.join('%.1f/%d/%.2f' % r for r in rows)))
