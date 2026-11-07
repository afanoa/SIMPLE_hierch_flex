#!/usr/bin/env python3
# continuous coordinate from a refine3D pass against K state maps: per particle the best dist to every state map (sparse
# prob_neigh part files dist_part_neigh_*.dat), centred per particle -> PCA over particles -> PC1 = coordinate; quantile bins
# of it = the next leaf states. Truth readout (--gt) for the benchmark only.
# usage: lrec_pc.py OUTDIR --tabs 'dir/dist_part_neigh_*.dat' --nstates K --nrows N --k KL [--prev coord.txt] [--gt gt.txt]
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--tabs', required=True); ap.add_argument('--nstates', type=int, required=True)
ap.add_argument('--nrows', type=int, required=True); ap.add_argument('--k', type=int, default=8); ap.add_argument('--prev'); ap.add_argument('--gt')
ap.add_argument('--npc', type=int, default=1); a = ap.parse_args()
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
pind = np.concatenate(P); D = np.concatenate(Ds).astype(float); full = np.isfinite(D).all(1)
print('tables: %d particles, %d with all %d states scored' % (len(pind), full.sum(), a.nstates))
pind, D = pind[full], D[full]; Dc = D - D.mean(1, keepdims=True); Dc -= Dc.mean(0)
U, S, Vt = np.linalg.svd(Dc, full_matrices=False); c = U[:, 0] * S[0]
if a.prev and os.path.exists(a.prev):
    pv = np.loadtxt(a.prev)[pind - 1]; ok = np.isfinite(pv)
    if np.corrcoef(c[ok], pv[ok])[0, 1] < 0: c = -c; Vt[0] = -Vt[0]
print('PCA variance fractions: %s; PC1 loadings: %s' % (' '.join('%.3f' % v for v in S[:5] ** 2 / (S ** 2).sum()), ' '.join('%.2f' % v for v in Vt[0])))
q = np.quantile(c, np.linspace(0, 1, a.k + 1)[1:-1]); b = np.searchsorted(q, c) + 1
os.makedirs(a.out, exist_ok=True); coord = np.full(a.nrows, np.nan); coord[pind - 1] = c; lab = np.zeros(a.nrows, int); lab[pind - 1] = b
np.savetxt(os.path.join(a.out, 'coord.txt'), coord, fmt='%.6g'); np.savetxt(os.path.join(a.out, 'bins.txt'), lab, fmt='%d')
for j in range(1, a.npc):
    pc = np.full(a.nrows, np.nan); pc[pind - 1] = U[:, j] * S[j]; np.savetxt(os.path.join(a.out, 'pc%d.txt' % (j + 1)), pc, fmt='%.6g')
if a.gt:
    gt = np.loadtxt(a.gt, comments='#', usecols=(1,), dtype=int)[pind - 1]; sw = (gt >= 200) & (gt <= 230); x = gt - 200.
    def eta2(lab_):
        W = sum(((x[sw & (lab_ == v)] - x[sw & (lab_ == v)].mean()) ** 2).sum() for v in np.unique(lab_[sw])); return 1 - W / ((x[sw] - x[sw].mean()) ** 2).sum()
    r = np.corrcoef(c[sw], x[sw])[0, 1]; s = 1 if r > 0 else -1; o = np.argsort(s * c); n = len(o)
    print('TRUTH: swing share %.3f, corr(PC1, swing) %.3f; K%d bins eta2 %.3f, bin means %s; bottom/top 2%% %.1f/%.1f, 5%% %.1f/%.1f'
          % (sw.mean(), abs(r), a.k, eta2(b), ' '.join('%.1f' % x[sw & (b == v)].mean() for v in (range(1, a.k + 1) if s > 0 else range(a.k, 0, -1))),
             x[o[:n // 50]].mean(), x[o[-n // 50:]].mean(), x[o[:n // 20]].mean(), x[o[-n // 20:]].mean()))
    for j in range(1, min(4, Dc.shape[1])): print('TRUTH: corr(PC%d, swing) %.3f' % (j + 1, abs(np.corrcoef(U[sw, j], x[sw])[0, 1])))
