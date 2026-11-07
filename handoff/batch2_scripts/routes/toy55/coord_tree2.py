#!/usr/bin/env python3
"""Stage 4 of the general pipeline, every branch at once: one cross-half scoring pass of all particles against all final-state maps
(refine3D prob tables, best dist per particle and state) -> for every node of the final-state tree (state_merge.py merge_tree.txt over
the final-state maps) holding >= MINP particles: the node's particles x the node's states' score columns, centred per particle and per
column -> PCA -> PC1 = the node's coordinate (PC1 variance share reported; a single final state is split on its own leaves' / a random
start's scores only when given) -> even-width sampling (2nd..98th percentile, end bins open, bins below MINB folded inward).
With GT + CODES: per node corr(PC1, truth coordinate), even-bin means, bottom / top 2 %.
usage: coord_tree.py OUT --tabs 'dir/dist_part_neigh_*.dat' --nstates K --units FINAL_UNITS --order 1,2,..(unit id of map 1..K)
       --tree merge_tree.txt [--minp 10000 --even 8 --minb 1000 --gt GT --codes swing:200:1]"""
import os, glob, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--tabs', required=True); ap.add_argument('--nstates', type=int, required=True)
ap.add_argument('--units', required=True); ap.add_argument('--order', required=True); ap.add_argument('--tree', required=True); ap.add_argument('--minp', type=int, default=10000)
ap.add_argument('--even', type=int, default=8); ap.add_argument('--dip', type=float, default=0.6); ap.add_argument('--pc1', type=float, default=0.65); ap.add_argument('--cum', type=float, default=0.70); ap.add_argument('--maxax', type=int, default=3); ap.add_argument('--grid', default='8,4,3'); ap.add_argument('--minb', type=int, default=1000); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default=''); ap.add_argument('--names', default='', help='file: unit id -> original state id (for --arm)'); ap.add_argument('--arm', default='', help='arm_states.txt of the original states: closed-arm readout per node')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
units = np.loadtxt(a.units, dtype=int); nrows = len(units); order = [int(x) for x in a.order.split(',')]   # map k (1-based) = unit order[k-1]
CAND = np.dtype([('iref', '<i4'), ('inpl', '<i4'), ('dist', '<f4'), ('x', '<f4'), ('y', '<f4'), ('has_sh', '<i4')]); P, Ds = [], []
for f in sorted(glob.glob(a.tabs)):   # the prob_neigh part tables, as lrec_pc.py
    sz = os.path.getsize(f); raw = open(f, 'rb').read(); nrefs, nptcls, nnz, nchunks = np.frombuffer(raw, '<i8', 4); nspace = nrefs // a.nstates
    pinds = np.frombuffer(raw, '<i4', nptcls, 32); pos = sz - (28 * nnz + 4 * nchunks); D = np.full((nptcls, a.nstates), np.inf, np.float32)
    for c in range(nchunks):
        n = int(np.frombuffer(raw, '<i4', 1, pos)[0]); pos += 4; idx = np.frombuffer(raw, '<i4', n, pos) - 1; pos += 4 * n
        cand = np.frombuffer(raw, CAND, n, pos); pos += 24 * n; st = (cand['iref'] - 1) // nspace
        ok = (cand['inpl'] > 0) & (st >= 0) & (st < a.nstates); np.minimum.at(D, (idx[ok], st[ok]), cand['dist'][ok])
    P.append(pinds.copy()); Ds.append(D)
pind = np.concatenate(P); D = np.concatenate(Ds).astype(float); print('score table: %d particles x %d states' % D.shape)
row = np.full(nrows + 1, -1); row[pind] = np.arange(len(pind))
hdr = open(a.tree).readline(); names = [int(x) for x in hdr.split('are the leaves ')[1].split(')')[0].split()]; N = len(names)
mem = {i: [names[i]] for i in range(N)}; S = np.atleast_2d(np.loadtxt(a.tree, comments='#')); nodes = []
for st in S: p, q, r = int(st[1]), int(st[2]), int(st[3]); mem[r] = mem[p] + mem[q]; nodes.append((r, mem[r]))
armv = {}; nm_ = {}
if a.names: nm_ = {int(l.split()[0]): int(l.split()[1]) for l in open(a.names) if l.strip() and not l.startswith('#')}
if a.arm:
    import re
    for l in open(a.arm):
        if l.strip() and not l.startswith('#'): w_ = l.split(); armv[int(re.findall(r'(\d+)\.mrc$', w_[0])[0])] = float(w_[1])
truth = {}
if a.gt:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(nrows, -1); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
    for cd in [w for w in a.codes.split(',') if w]:
        p = cd.split(':')
        if p[0] == 'grid': truth['grid_k'] = np.where(conf >= int(p[1]), ((conf - int(p[1])) // 10) * float(p[2]), np.nan); truth['grid_j'] = np.where(conf >= int(p[1]), ((conf - int(p[1])) % 10) * float(p[3]), np.nan)
        else:
            hi = min([int(q.split(':')[1]) for q in a.codes.split(',') if q and not q.startswith('grid:') and int(q.split(':')[1]) > int(p[1])] + [10 ** 9]); truth[p[0]] = np.where((conf >= int(p[1])) & (conf < hi), (conf - int(p[1])) * float(p[2]), np.nan)
def dip_ratio(c):   # deepest interior valley of the (3-bin smoothed, 24-bin) coordinate histogram relative to the smaller flanking peak
    lo_, hi_ = np.percentile(c, [1, 99]); h = np.histogram(np.clip(c, lo_, hi_), 24, (lo_, hi_))[0].astype(float); h = np.convolve(h, np.ones(3) / 3, 'same'); best = 1.
    for i in range(1, 23): best = min(best, h[i] / max(min(h[:i].max(), h[i + 1:].max()), 1e-9))
    return best
cont = {}; parent = {}
for st in S: parent[int(st[1])] = int(st[3]); parent[int(st[2])] = int(st[3])
for r, ms in nodes:
    cols = [order.index(m) for m in ms]; rws = np.where(np.isin(units, ms))[0]; rws = rws[row[rws + 1] >= 0]
    if len(rws) < a.minp: continue
    X = D[row[rws + 1]][:, cols]; ok = np.isfinite(X).all(1); rws, X = rws[ok], X[ok]; Xc = X - X.mean(1, keepdims=True); Xc -= Xc.mean(0)
    U, Sv, Vt = np.linalg.svd(Xc, full_matrices=False); c = U[:, 0] * Sv[0]; sh = Sv ** 2 / (Sv ** 2).sum(); share = sh[0]
    cs = np.cumsum(sh); k = int(np.searchsorted(cs, a.cum) + 1)   # fewest axes carrying CUM of the node's score variance
    if len(ms) == 2: k = 1
    Zc = U[:, :max(k, 1)] * Sv[:max(k, 1)]; G_ = [int(g) for g in a.grid.split(',')]; nb = G_[min(k, len(G_)) - 1] if k <= a.maxax else 0
    cells = np.zeros(len(c), int)
    if nb:
        idx = []
        for j in range(k):
            lo, hi = np.percentile(Zc[:, j], [2, 98]); ed = np.linspace(lo, hi, nb + 1); idx.append(np.clip(np.searchsorted(ed, Zc[:, j], side='right') - 1, 0, nb - 1))
        cells = np.ravel_multi_index(idx, [nb] * k)
        sz = np.bincount(cells, minlength=nb ** k); big = np.where(sz >= a.minb)[0]
        if len(big):
            cen = np.array([Zc[cells == q].mean(0) if sz[q] else np.zeros(k) for q in range(nb ** k)])
            for q in np.where((sz > 0) & (sz < a.minb))[0]: cells[cells == q] = big[np.argmin(((cen[big] - cen[q]) ** 2).sum(1))]   # small cells fold into the nearest kept cell
        _, cells = np.unique(cells, return_inverse=True)
    segs = [list(np.where(cells == q)[0]) for q in range(cells.max() + 1)]
    lab = np.zeros(nrows, int); co = np.full(nrows, np.nan); co[rws] = c
    for q_, g in enumerate(segs, 1): lab[rws[g]] = q_
    d = '%s/node%03d' % (a.out, r); os.makedirs(d, exist_ok=True); np.savetxt(d + '/coord.txt', co, fmt='%.6g'); np.savetxt(d + '/even_bins.txt', lab, fmt='%d')
    for j in range(1, Zc.shape[1]):
        cj = np.full(nrows, np.nan); cj[rws] = Zc[:, j]; np.savetxt(d + '/coord_axis%d.txt' % (j + 1), cj, fmt='%.6g')
    drs = [dip_ratio(Zc[:, j]) for j in range(min(k, Zc.shape[1]))]; dr = min(drs); cont[r] = bool(nb) and dr >= a.dip   # <= MAXAX axes carry CUM and none has a deep valley
    line = 'node %d: states %s, %d particles; axis shares %s -> %d axes; valley ratio %.2f (%s); %d samples (%s)' % (r, '+'.join(map(str, ms)), len(rws), ' '.join('%.2f' % v for v in sh[:4]), k, dr, 'samplable' if cont[r] else ('gapped' if dr < a.dip else 'too many axes'), len(segs), ' '.join(str(len(g)) for g in segs))
    if armv:
        av = [(armv.get(nm_.get(m, m), np.nan), (units == m).sum(), nm_.get(m, m)) for m in ms]; av = [t for t in av if np.isfinite(t[0])]
        if av: line += '\n      closed-arm of the member states (original ids): weighted mean %.2f, range %.2f..%.2f; %s' % (sum(x * n for x, n, _ in av) / sum(n for _, n, _ in av), min(x for x, _, _ in av), max(x for x, _, _ in av), ' '.join('%d:%.2f' % (o, x) for x, _, o in sorted(av, key=lambda t: t[0])))
    for nm, x in truth.items():
        xs = x[rws]; sw = np.isfinite(xs)
        if sw.sum() < 0.3 * len(rws) or sw.sum() < 200: continue
        rr = np.corrcoef(c[sw], xs[sw])[0, 1]; sg = 1 if rr > 0 else -1; o = np.argsort(sg * c[sw]); n = len(o); lb = lab[rws][sw]
        em = [xs[sw][lb == v].mean() for v in (range(1, len(segs) + 1) if sg > 0 else range(len(segs), 0, -1)) if (lb == v).any()]
        line += '\n      %s (%.0f %% of the node): corr %.2f, even-bin means %s, bottom/top 2 %% %.1f/%.1f' % (nm, 100 * sw.mean(), abs(rr), ' '.join('%.1f' % v for v in em), xs[sw][o[:max(1, n // 50)]].mean(), xs[sw][o[-max(1, n // 50):]].mean())
    print(line)
# sampled branches: the maximal continuous nodes (continuous, parent gapped or too small)
sel = [r for r in cont if cont[r] and not cont.get(parent.get(r, -1), False)]
print('== sampled branches (maximal samplable nodes): %s' % (', '.join('node %d' % r for r in sorted(sel)) or 'none'))
open(os.path.join(a.out, 'sampled_nodes.txt'), 'w').write(''.join('%d\n' % r for r in sorted(sel)))
