#!/usr/bin/env python3
"""land_basins.py: BASIN TREE of an occupancy landscape along a path (energy-landscape disconnectivity graph): p(s) from the validated
reweighting (land_pscore lmarg likelihoods, pixel-split temperature, smoothing LAM), fitted on all particles and on each half separately.
Basins = local maxima of p; neighbouring basins join at the saddle between them; the join with the lowest barrier (saddle G minus the
shallower basin's G) is merged repeatedly until every remaining barrier exceeds THR kT; each merge is an internal node (its height = the
barrier), the surviving basins are the leaves. A split counts only if the same number of basins survives in BOTH halves with each
boundary within TOL path steps. Reports basins (position, mass, G), merge heights, and the half check, for each threshold.
Synthetic mode (--gt): the same on imposed shapes. usage: land_basins.py --scores NPZ [--gt GT] [--lam 0.01 --thr 0.3,0.5,0.75 --exclude F]"""
import os, sys, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--gt', default=''); ap.add_argument('--ids', default='200:230')
ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder'); ap.add_argument('--variant', default='lmarg'); ap.add_argument('--tau', default='split')
ap.add_argument('--lam', type=float, default=0.01); ap.add_argument('--thr', default='0.3,0.5,0.75'); ap.add_argument('--tol', type=int, default=3); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--exclude', default='', help='file of master-row flags (nonzero = drop), e.g. junk leaves'); a = ap.parse_args()
def log(*s): print(*s, flush=True)
Z = np.load(a.scores); rows = Z['rows']; K = Z['l0'].shape[1]; tg = np.arange(K, dtype=float); tau = float(np.median(Z['tau_split'])) if a.tau == 'split' else float(a.tau)
half = Z['eo'].astype(int) if ((Z['eo'] == 0).any() and (Z['eo'] == 1).any()) else Z['half'].astype(int)
keep = np.ones(len(rows), bool)
if a.exclude: drop = np.loadtxt(a.exclude).astype(int) != 0; keep = ~drop[rows]; log('excluded %d of %d particles (%s)' % ((~keep).sum(), len(rows), a.exclude))
def fit(idx):
    l = Z[a.variant][idx].astype(np.float64) / tau; l -= l.max(1, keepdims=True); L = np.exp(l); return LD.npmle(L) if a.lam == 0 else LD.fit_pen(L, a.lam)[1]
def basins(p, thr):
    G = -np.log(np.maximum(p, 1e-12)); mx = [i for i in range(K) if (i == 0 or G[i] <= G[i - 1]) and (i == K - 1 or G[i] <= G[i + 1])]
    B = [dict(lo=0, hi=K - 1, m=mx[0])]   # basins as [lo, hi] catchments between saddles, each with its minimum m
    for j in range(1, len(mx)):
        sd = mx[j - 1] + int(np.argmax(G[mx[j - 1]:mx[j] + 1])); B[-1]['hi'] = sd; B.append(dict(lo=sd + 1, hi=K - 1, m=mx[j]))
    merges = []
    while len(B) > 1:
        bar = [G[B[k]['hi']] - max(G[B[k]['m']], G[B[k + 1]['m']]) for k in range(len(B) - 1)]; k = int(np.argmin(bar))
        if bar[k] >= thr: break
        merges.append((round(float(bar[k]), 2), B[k]['m'], B[k + 1]['m'])); m = B[k]['m'] if G[B[k]['m']] <= G[B[k + 1]['m']] else B[k + 1]['m']
        B[k:k + 2] = [dict(lo=B[k]['lo'], hi=B[k + 1]['hi'], m=m)]
    for b in B: b['mass'] = float(p[b['lo']:b['hi'] + 1].sum()); b['G'] = float(G[b['m']] - G.min())
    bars = [float(G[B[k]['hi']] - max(G[B[k]['m']], G[B[k + 1]['m']])) for k in range(len(B) - 1)]
    return B, bars, merges
def report(lab, idx, ptrue=None):
    pf = fit(idx); ph = [fit(idx[half[idx] == h]) for h in (0, 1)]
    for thr in [float(x) for x in a.thr.split(',')]:
        Bf, bf, mf = basins(pf, thr); Bh = [basins(q, thr) for q in ph]
        same = all(len(b[0]) == len(Bf) for b in Bh) and all(all(abs(x['hi'] - y['hi']) <= a.tol for x, y in zip(b[0][:-1], Bf[:-1])) for b in Bh)
        s_ = ' | '.join('basin s %.2f (%.2f-%.2f) mass %.3f G %.2f' % (b['m'] / (K - 1.), b['lo'] / (K - 1.), b['hi'] / (K - 1.), b['mass'], b['G']) for b in Bf)
        tr = ''
        if ptrue is not None: Bt, bt, _ = basins(ptrue, thr); tr = ' || TRUTH %d basin(s)%s' % (len(Bt), (' barriers ' + ' '.join('%.2f' % x for x in bt)) if bt else '')
        log('RESULT %-10s thr %.2f: %d basin(s) [halves %d / %d, %s]%s | %s | merged below thr: %s%s' % (lab, thr, len(Bf), len(Bh[0][0]), len(Bh[1][0]), 'REPRODUCED' if same else 'not reproduced',
            (' barriers ' + ' '.join('%.2f' % x for x in bf)) if bf else '', s_, ' '.join('%.2f@%d-%d' % m for m in mf) if mf else '-', tr))
if a.gt:
    i0, i1 = [int(x) for x in a.ids.split(':')]; conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int); tt = conf[rows] - i0; ok = (tt >= 0) & (tt < K) & keep; Nt = np.bincount(tt[ok], minlength=K)
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
        if p is None: idx = cod
        else:
            sc = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(sc * p).astype(int); r = np.random.default_rng(a.seed); idx = np.concatenate([r.choice(cod[tt[cod] == k], nk[k], replace=False) for k in range(K)])
        report(s, idx, np.bincount(tt[idx], minlength=K) / len(idx))
else:
    report('real', np.where(keep)[0])
log('DONE')
