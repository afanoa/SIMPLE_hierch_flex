#!/usr/bin/env python3
"""land_twostate.py: particle-level FEW-STATE vs PATH test from land_pscore.py likelihoods. Mixtures over (a) the best TWO path maps at free
positions, (b) the best THREE, (c) every path map (NPMLE); positions and weights fitted on one half and scored on the other (both ways); the
held-out gain of (c) over (a) and (b) per particle (nats), with a particle bootstrap. Particles that sit between the ends need the intermediate
maps: a gain above the two-state null means intermediates are populated. Synthetic (--gt): the same on imposed shapes, which calibrates
the gain (dw4 is nearly two-state, flat is a continuum). usage: land_twostate.py --scores NPZ [--variant lmarg --tau split] [--gt GT]"""
import os, sys, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--variant', default='lmarg'); ap.add_argument('--tau', default='split')
ap.add_argument('--gt', default=''); ap.add_argument('--ids', default='200:230'); ap.add_argument('--shapes', default='dw4,dw2,flat,broad'); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--nboot', type=int, default=200); ap.add_argument('--exclude', default=''); a = ap.parse_args()
def log(*s): print(*s, flush=True)
Z = np.load(a.scores); K = Z['l0'].shape[1]; tau = float(np.median(Z['tau_split'])) if a.tau == 'split' else float(a.tau); rows = Z['rows']; tg = np.arange(K, dtype=float)
half = Z['eo'].astype(int) if ((Z['eo'] == 0).any() and (Z['eo'] == 1).any()) else Z['half'].astype(int)
if half.min() == half.max(): half = np.random.default_rng(a.seed).integers(0, 2, len(rows))
def em_w(L, nit=200):
    w = np.full(L.shape[1], 1. / L.shape[1])
    for _ in range(nit): d = L @ w; w = w * (L.T @ (1. / d)) / len(L)
    return w
def best_subset(Ltr, k):   # greedy-exhaustive: best pair over all pairs; best triple = best pair + the best third map
    if k == 2:
        best = None
        for i in range(K):
            for j in range(i + 1, K):
                w = em_w(Ltr[:, [i, j]], 100); v = np.log(Ltr[:, [i, j]] @ w).mean()
                if best is None or v > best[0]: best = (v, [i, j])
        return best[1]
    pr = best_subset(Ltr, 2); best = None
    for m in range(K):
        if m in pr: continue
        c = sorted(pr + [m]); w = em_w(Ltr[:, c], 100); v = np.log(Ltr[:, c] @ w).mean()
        if best is None or v > best[0]: best = (v, c)
    return best[1]
def heldout_gain(idx):
    l = Z[a.variant][idx].astype(np.float64) / tau; l -= l.max(1, keepdims=True); L = np.exp(l); h = half[idx]; out = {}; picks = {}
    for name in ('two', 'three', 'path'):
        ll = np.zeros(len(idx)); picks[name] = []
        for g in (0, 1):
            Ltr, Lte = L[h == g], L[h != g]
            comp = list(range(K)) if name == 'path' else best_subset(Ltr, 2 if name == 'two' else 3); picks[name].append(comp if name != 'path' else 'all')
            p = LD.npmle(Ltr[:, comp]) if name == 'path' else em_w(Ltr[:, comp]); ll[h != g] = np.log(np.maximum(Lte[:, comp] @ p, 1e-300))
        out[name] = ll
    d = out['path'] - out['two']; d3 = out['path'] - out['three']; rng = np.random.default_rng(a.seed); bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(a.nboot)]
    bs3 = [d3[rng.integers(0, len(d3), len(d3))].mean() for _ in range(a.nboot)]
    log('   picks: two %s, three %s' % (picks['two'], picks['three']))
    return d.mean(), np.percentile(bs, [2.5, 97.5]), d3.mean(), np.percentile(bs3, [2.5, 97.5])
if a.gt:
    i0, i1 = [int(x) for x in a.ids.split(':')]; conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int); tt = conf[rows] - i0; ok = (tt >= 0) & (tt < K); Nt = np.bincount(tt[ok], minlength=K)
    def target(name):
        if name == 'flat': g = np.zeros(K)
        elif name.startswith('dw'): hh = float(name[2:]); g = hh * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
        elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
        else: return None
        p = np.exp(-g); return p / p.sum()
    for s in a.shapes.split(','):
        cod = np.where(ok)[0]; p = target(s)
        if p is None: idx = cod
        else:
            sc = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(sc * p).astype(int); r = np.random.default_rng(a.seed); idx = np.concatenate([r.choice(cod[tt[cod] == k], nk[k], replace=False) for k in range(K)])
        m, ci, m3, ci3 = heldout_gain(idx)
        log('RESULT %-6s n %6d: held-out gain of the path over the best two states %.4f nats/ptcl [%.4f %.4f], over the best three %.4f [%.4f %.4f]' % (s, len(idx), m, ci[0], ci[1], m3, ci3[0], ci3[1]))
else:
    sel = np.arange(len(rows)) if not a.exclude else np.where(np.loadtxt(a.exclude).astype(int)[rows] == 0)[0]
    m, ci, m3, ci3 = heldout_gain(sel)
    log('RESULT real %s (%d particles, %s tau %.3f): held-out gain of the path over the best two states %.4f nats/ptcl [%.4f %.4f], over the best three %.4f [%.4f %.4f]' % (
        a.scores, len(rows), a.variant, tau, m, ci[0], ci[1], m3, ci3[0], ci3[1]))
