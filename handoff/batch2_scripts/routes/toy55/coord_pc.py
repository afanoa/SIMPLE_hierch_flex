#!/usr/bin/env python3
# coord_pc.py = lrec_pc.py (realistic-SNR session) + --codes truth readout for any panel set + even-width bins of the coordinate
# (EVEN: 2nd..98th percentile, equal width, end bins open, bins below MINB particles folded inward = the sampled trajectory).
# continuous coordinate from a refine3D pass against K state maps: per particle the best dist to every state map (sparse
# prob_neigh part files dist_part_neigh_*.dat), centred per particle -> PCA over particles -> PC1 = coordinate; quantile bins
# of it = the next leaf states. Truth readout (--gt) for the benchmark only.
# usage: lrec_pc.py OUTDIR --tabs 'dir/dist_part_neigh_*.dat' --nstates K --nrows N --k KL [--prev coord.txt] [--gt gt.txt]
import argparse, glob, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--tabs', required=True); ap.add_argument('--nstates', type=int, required=True)
ap.add_argument('--nrows', type=int, required=True); ap.add_argument('--k', type=int, default=8); ap.add_argument('--prev'); ap.add_argument('--gt')
ap.add_argument('--npc', type=int, default=1); ap.add_argument('--codes', default='swing:200:1'); ap.add_argument('--even', type=int, default=8); ap.add_argument('--minb', type=int, default=1000); a = ap.parse_args()
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
# even-width sampling of the coordinate: equal width between the 2nd and 98th percentile, end bins open, small bins folded inward
lo, hi = np.percentile(c, [2, 98]); ed = np.linspace(lo, hi, a.even + 1); eb = np.clip(np.searchsorted(ed, c, side='right') - 1, 0, a.even - 1)
segs = [list(np.where(eb == k)[0]) for k in range(a.even)]; segs = [g for g in segs if g]
for e in (0, -1):
    while len(segs) > 2 and len(segs[e]) < a.minb: k = 1 if e == 0 else -2; segs[k] = segs[k] + segs[e]; segs.pop(e)
i = 1
while i < len(segs) - 1:
    if len(segs[i]) < a.minb: j = i - 1 if len(segs[i - 1]) <= len(segs[i + 1]) else i + 1; segs[j] = segs[j] + segs[i]; segs.pop(i)
    else: i += 1
ev_ = np.zeros(len(c), int)
for k, g in enumerate(segs, 1): ev_[g] = k
elab = np.zeros(a.nrows, int); elab[pind - 1] = ev_; np.savetxt(os.path.join(a.out, 'even_bins.txt'), elab, fmt='%d')
print('even sampling: %d bins, particles %s' % (len(segs), ' '.join(str(len(g)) for g in segs)))
if a.gt:
    gt = np.loadtxt(a.gt, comments='#', usecols=(1,), dtype=int)[pind - 1]
    for cd in [w for w in a.codes.split(',') if w]:
        p = cd.split(':')
        if p[0] == 'grid': tv = {'grid_k': np.where(gt >= int(p[1]), ((gt - int(p[1])) // 10) * float(p[2]), np.nan), 'grid_j': np.where(gt >= int(p[1]), ((gt - int(p[1])) % 10) * float(p[3]), np.nan)}
        else:
            hi_ = min([int(q.split(':')[1]) for q in a.codes.split(',') if q and not q.startswith('grid:') and int(q.split(':')[1]) > int(p[1])] + [10 ** 9])
            tv = {p[0]: np.where((gt >= int(p[1])) & (gt < hi_), (gt - int(p[1])) * float(p[2]), np.nan)}
        for nm, x in tv.items():
            sw = np.isfinite(x)
            if sw.sum() < 50: print('TRUTH %s: %d particles of it in this node' % (nm, sw.sum())); continue
            r = np.corrcoef(c[sw], x[sw])[0, 1]; sg = 1 if r > 0 else -1; o = np.argsort(sg * c[sw]); n = len(o); xs = x[sw]
            def eta2(lab_):
                W = sum(((xs[lab_[sw] == v] - xs[lab_[sw] == v].mean()) ** 2).sum() for v in np.unique(lab_[sw])); return 1 - W / max(((xs - xs.mean()) ** 2).sum(), 1e-12)
            em = [xs[ev_[sw] == v].mean() for v in (range(1, len(segs) + 1) if sg > 0 else range(len(segs), 0, -1)) if (ev_[sw] == v).any()]
            print('TRUTH %s: share %.3f of the node, corr(PC1, %s) %.3f; quantile bins eta2 %.3f; even bins eta2 %.3f, means %s; bottom/top 2%% %.1f/%.1f'
                  % (nm, sw.mean(), nm, abs(r), eta2(b), eta2(ev_), ' '.join('%.1f' % v for v in em), xs[o[:max(1, n // 50)]].mean(), xs[o[-max(1, n // 50):]].mean()))
