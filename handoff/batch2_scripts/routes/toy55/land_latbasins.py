#!/usr/bin/env python3
"""land_latbasins.py: BASINS of the particle density in a flex latent space, with no path and no K. kNN log-density at every particle
(ln f = -D ln r_k), then topological mode analysis (ToMATo, Chazal et al. 2013) on the symmetric kNN graph: points enter in order of
decreasing density and join their densest earlier neighbour's cluster; when a point links two clusters, the lower-peaked one merges if
its prominence (its peak ln f minus the linking point's ln f, in kT) is below THR. Surviving clusters are the basins (every particle is
assigned to one); the prominences of all clusters form the 0-dim persistence diagram. A split counts only if each half (eo) on its own
gives the same number of basins and its partition matches the full one (ARI >= 0.8 on the half's particles).
Modes:
  --gt G --codes 200:230 --coords C : imposed shapes along a known coordinate (swing30 t) by subsampling, as land_basins.py
  --gt G --discrete --coords C      : basins against discrete truth groups (contingency), e.g. simsites combinations
  --coords A,B,C --names a,b,c      : basins per run and per half, then cross-run ARI and matched occupancies; --labout writes, per
                                      run and threshold, one basin per master row (0 = not in the node) for land_rec.sbatch
usage: land_latbasins.py --coords F[,F2..] [--dims 10 --k 200 --kg 20 --thr 0.3,0.5,0.75,1.0,1.5 --gt G --codes I0:I1 | --discrete]"""
import argparse, numpy as np
from scipy.spatial import cKDTree
ap = argparse.ArgumentParser(); ap.add_argument('--coords', required=True); ap.add_argument('--names', default=''); ap.add_argument('--dims', default='10')
ap.add_argument('--k', default='200'); ap.add_argument('--kg', type=int, default=20); ap.add_argument('--thr', default='0.5,0.75,1.0,1.5')
ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='200:230'); ap.add_argument('--discrete', action='store_true')
ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder'); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--exclude', default='', help='file of master-row flags (nonzero = drop)'); ap.add_argument('--labout', default='')
ap.add_argument('--labthr', default='0.5,1.0'); ap.add_argument('--selftest', action='store_true'); ap.add_argument('--nm', type=int, default=0, help='master rows for --labout files'); a = ap.parse_args()
def log(*s): print(*s, flush=True)
THR = [float(x) for x in a.thr.split(',')]; dims = [int(x) for x in a.dims.split(',')]; ks = [int(x) for x in a.k.split(',')]

def load(path):
    X = np.loadtxt(path, comments='#'); return X[:, 0].astype(int) - 1, X[:, 1].astype(int), X[:, 5:]
def knn(z, k, kg):
    d, i = cKDTree(z).query(z, k=max(k, kg) + 1, workers=-1)
    return -z.shape[1] * np.log(np.maximum(d[:, k], 1e-12)), i[:, 1:kg + 1]
def tomato(lnf, nbr, tau):
    """ToMATo clustering; returns root labels (0..nb-1 by decreasing peak), prominences of every cluster that merged (tau = inf
    gives the full diagram), and the number of graph components."""
    n = len(lnf); order = np.argsort(-lnf, kind='stable'); rank = np.empty(n, np.int64); rank[order] = np.arange(n)
    a_ = np.repeat(np.arange(n), nbr.shape[1]); b_ = nbr.ravel(); A = np.concatenate([a_, b_]); B = np.concatenate([b_, a_])
    m = rank[B] < rank[A]; key = np.unique(A[m].astype(np.int64) * n + B[m]); A = (key // n).astype(np.int64); B = (key % n).astype(np.int64)
    start = np.searchsorted(A, np.arange(n + 1)); Bl = B.tolist(); st = start.tolist(); f = lnf.tolist(); parent = list(range(n)); pers = []
    def find(x):
        r = x
        while parent[r] != r: r = parent[r]
        while parent[x] != r: parent[x], x = r, parent[x]
        return r
    for i in order.tolist():
        nb = Bl[st[i]:st[i + 1]]
        if not nb: continue
        g = max(nb, key=lambda j: f[j]); parent[i] = find(g)
        for j in nb:
            rj = find(j); ri = find(i)
            if rj == ri: continue
            lo, hi = (rj, ri) if f[rj] < f[ri] else (ri, rj); prom = f[lo] - f[i]
            if prom < tau: parent[lo] = hi; pers.append(prom)
    roots = np.array([find(x) for x in range(n)]); ur = np.unique(roots); ur = ur[np.argsort(-lnf[ur])]
    lab = np.searchsorted(np.sort(ur), roots); remap = np.empty(len(ur), int); remap[np.searchsorted(np.sort(ur), ur)] = np.arange(len(ur))
    return remap[lab], np.sort(np.array(pers))[::-1], len(ur) if tau == np.inf else None
def ari(x, y):
    from math import comb
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); C = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(C, (xi, yi), 1)
    s = lambda v: (v * (v - 1) / 2.).sum(); n = len(x); t = s(C); ta = s(C.sum(1)); tb = s(C.sum(0)); e = ta * tb / (n * (n - 1) / 2.)
    return 1. if ta == tb == t else (t - e) / (0.5 * (ta + tb) - e)
def basins(z, k, kg):
    """full persistence diagram plus a labelling per threshold"""
    lnf, nbr = knn(z, k, kg); _, pinf, ncomp = tomato(lnf, nbr, np.inf)
    return lnf, nbr, pinf, ncomp, {t: tomato(lnf, nbr, t)[0] for t in THR}
def analyse(lab, z, eo, k, kg, ptrue=None, groups=None):
    """basins of z (all) and of each half; returns per-threshold labels of the full set"""
    lnf, nbr, pinf, ncomp, L = basins(z, k, kg); H = [basins(z[eo == h], max(10, k // 2), kg) for h in (0, 1)]
    gp = pinf[:9] - pinf[1:10] if len(pinf) > 1 else np.zeros(1); j = int(np.argmax(gp))
    log('  %-9s n %6d D %2d k %4d: components %d, top prominences %s (largest gap below #%d -> %d basins) | halves %s / %s' % (lab, len(z), z.shape[1], k, ncomp,
        ' '.join('%.2f' % x for x in pinf[:6]), j + 1, j + 1 + ncomp, ' '.join('%.2f' % x for x in H[0][2][:4]), ' '.join('%.2f' % x for x in H[1][2][:4])))
    for t in THR:
        lf = L[t]; nf = int(lf.max()) + 1; nh = [int(H[h][4][t].max()) + 1 for h in (0, 1)]
        ah = [ari(lf[eo == h], H[h][4][t]) for h in (0, 1)]; rep = (nh[0] == nf == nh[1]) and min(ah) >= 0.8
        ms = np.bincount(lf, minlength=nf) / len(lf); tr = ''
        if ptrue is not None: tr = ' || truth (1-D) %d basin(s)' % ptrue(t)
        if groups is not None:
            G, gl = groups; tr = ' || ARI vs truth groups %.3f; basin x group counts: %s' % (ari(lf, G), ' ; '.join('b%d: ' % b + ' '.join('%s=%d' % (gl[g], ((lf == b) & (G == g)).sum())
                for g in np.unique(G) if ((lf == b) & (G == g)).sum() > 0.01 * (lf == b).sum()) for b in range(min(nf, 8))))
        log('RESULT %-9s D %2d k %4d thr %.2f: %d basin(s) [halves %d / %d, ARI %.2f / %.2f, %s] masses %s%s' % (lab, z.shape[1], k, t, nf, nh[0], nh[1], ah[0], ah[1],
            'REPRODUCED' if rep else 'not reproduced', ' '.join('%.3f' % x for x in ms[:8]), tr))
    return L, pinf

def bas1d(p, thr):   # the 1-D basin rule of land_basins.py on a histogram (truth count)
    K = len(p); G = -np.log(np.maximum(p, 1e-12)); mx = [i for i in range(K) if (i == 0 or G[i] <= G[i - 1]) and (i == K - 1 or G[i] <= G[i + 1])]
    B = [dict(lo=0, hi=K - 1, m=mx[0])]
    for j in range(1, len(mx)):
        sd = mx[j - 1] + int(np.argmax(G[mx[j - 1]:mx[j] + 1])); B[-1]['hi'] = sd; B.append(dict(lo=sd + 1, hi=K - 1, m=mx[j]))
    while len(B) > 1:
        bar = [G[B[k]['hi']] - max(G[B[k]['m']], G[B[k + 1]['m']]) for k in range(len(B) - 1)]; k = int(np.argmin(bar))
        if bar[k] >= thr: break
        m = B[k]['m'] if G[B[k]['m']] <= G[B[k + 1]['m']] else B[k + 1]['m']; B[k:k + 2] = [dict(lo=B[k]['lo'], hi=B[k + 1]['hi'], m=m)]
    return len(B)

if a.selftest:   # one blob / two blobs 4 sd apart / a noisy segment, D 10, n 40000
    r = np.random.default_rng(1); n = 40000; eo_ = r.integers(0, 2, n)
    for lab_, z_ in [('blob', r.normal(size=(n, 10))), ('twoblobs', r.normal(size=(n, 10)) + np.outer(np.where(r.random(n) < 0.5, -2., 2.), np.eye(10)[0])),
                     ('segment', r.normal(size=(n, 10)) + np.outer(r.uniform(-4, 4, n), np.eye(10)[0]))]:
        for k in ks: analyse(lab_, z_, eo_, k, a.kg)
    raise SystemExit
paths = a.coords.split(','); names = a.names.split(',') if a.names else ['run%d' % i for i in range(len(paths))]
if a.gt:
    rows, eo, Z = load(paths[0]); gt = np.loadtxt(a.gt, comments='#'); conf = gt[:, 1].astype(int)[rows]
    keep = np.ones(len(rows), bool)
    if a.exclude: keep &= np.loadtxt(a.exclude).astype(int)[rows] == 0
    log('%s: %d particles, latent dims %d; per-dim sd %s' % (paths[0], len(rows), Z.shape[1], ' '.join('%.2f' % x for x in Z.std(0))))
    if a.discrete:
        ok = keep & (conf > 0); log('discrete truth: %d non-junk of %d; groups %s' % (ok.sum(), len(ok), ' '.join('%d:%d' % (g, (conf[ok] == g).sum()) for g in np.unique(conf[ok]))))
        for D in dims:
            for k in ks: analyse('nonjunk', Z[ok, :D], eo[ok], k, a.kg, groups=(conf[ok], {g: str(g) for g in np.unique(conf[ok])}))
            for k in ks: analyse('all', Z[keep, :D], eo[keep], k, a.kg, groups=(conf[keep], {g: str(g) for g in np.unique(conf[keep])}))
    else:
        i0, i1 = [int(x) for x in a.codes.split(':')]; K = i1 - i0 + 1; tt = conf - i0; ok = keep & (tt >= 0) & (tt < K); Nt = np.bincount(tt[ok], minlength=K); tg = np.arange(K, dtype=float)
        log('coordinate truth: %d of %d particles carry codes %d-%d; per-bin counts %d-%d' % (ok.sum(), len(ok), i0, i1, Nt.min(), Nt.max()))
        def target(name):
            if name == 'flat': g = np.zeros(K)
            elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
            elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
            elif name == 'shoulder':
                m_ = np.exp(-0.5 * ((tg - 20.) / 4.) ** 2); p = 0.95 * m_ / m_.sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
            else: return None
            p = np.exp(-g); return p / p.sum()
        for s in a.shapes.split(','):
            cod = np.where(ok)[0]; p = target(s)
            if s == 'node': idx = np.where(keep)[0]
            elif p is None: idx = cod
            else:
                sc = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(sc * p).astype(int); r = np.random.default_rng(a.seed)
                idx = np.concatenate([r.choice(cod[tt[cod] == kk], nk[kk], replace=False) for kk in range(K)])
            ph = np.bincount(tt[idx][tt[idx] >= 0][tt[idx][tt[idx] >= 0] < K], minlength=K) / max(1, len(idx))
            for D in dims:
                for k in ks: analyse(s, Z[idx, :D], eo[idx], k, a.kg, ptrue=(lambda t, ph=ph: bas1d(ph, t)))
else:
    R = {}
    for nm_, pth in zip(names, paths):
        rows, eo, Z = load(pth); keep = np.ones(len(rows), bool)
        if a.exclude: keep &= np.loadtxt(a.exclude).astype(int)[rows] == 0
        log('%s: %s, %d particles (%d kept), per-dim sd %s' % (nm_, pth, len(rows), keep.sum(), ' '.join('%.2f' % x for x in Z.std(0))))
        R[nm_] = dict(rows=rows[keep], eo=eo[keep], Z=Z[keep])
    common = sorted(set.intersection(*[set(r['rows'].tolist()) for r in R.values()])); log('particles common to all runs: %d' % len(common))
    for D in dims:
        for k in ks:
            LB = {}
            for nm_, r in R.items():
                L, _ = analyse(nm_, r['Z'][:, :D], r['eo'], k, a.kg); LB[nm_] = {t: dict(zip(r['rows'].tolist(), L[t].tolist())) for t in THR}
                if a.labout and a.nm:
                    for t in [float(x) for x in a.labthr.split(',')]:
                        if t not in L: continue
                        lab = np.zeros(a.nm, int); lab[r['rows']] = L[t] + 1; np.savetxt('%s_%s_D%d_k%d_thr%.2f.txt' % (a.labout, nm_, D, k, t), lab, fmt='%d')
            for t in THR:
                nn = list(LB); M = {x: np.array([LB[x][t][c] for c in common]) for x in nn}
                pr = ' '.join('%s-%s %.3f' % (x, y, ari(M[x], M[y])) for i, x in enumerate(nn) for y in nn[i + 1:])
                occ = ' | '.join('%s %s' % (x, ' '.join('%.3f' % v for v in (np.bincount(M[x]) / len(common))[:8])) for x in nn)
                # matched occupancies: each basin of the first run against its best-overlap basin in every other run
                ref = nn[0]; mt = []
                for b in range(M[ref].max() + 1):
                    sel = M[ref] == b; row_ = ['%s b%d %.3f' % (ref, b, sel.mean())]
                    for y in nn[1:]:
                        cnt = np.bincount(M[y][sel], minlength=M[y].max() + 1); bb = int(np.argmax(cnt)); jac = cnt[bb] / float((sel | (M[y] == bb)).sum())
                        row_.append('%s b%d %.3f (Jaccard %.2f)' % (y, bb, (M[y] == bb).mean(), jac))
                    mt.append(', '.join(row_))
                log('CROSS D %2d k %4d thr %.2f: ARI %s | occupancies %s | matched: %s' % (D, k, t, pr, occ, ' ;; '.join(mt[:6])))
log('DONE')
