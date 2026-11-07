#!/usr/bin/env python3
"""land_e0.py: ORACLE information ceiling of an occupancy landscape from static images. Every swing particle is scored against the TRUE
swing maps (lsamp_score2.py with the truth maps as leaves). Known landscapes p(t) are imposed by subsampling the particles by their true
swing t, and p(t) is recovered by (1) hard assignment (argmax), (2) a per-particle point coordinate (posterior mean under a flat prior,
histogrammed: what a latent embedding plus a density gives), (3) the NPMLE of the grid weights (EM over the per-particle likelihoods,
SQUAREM), with bootstrap bands. Scored: total variation to the truth, the low-t basin mass, the minor (t <= 3) mass, and the barrier
G(15) - min G(5), G(25) on 3-bin averages (G = -ln p, kT units).
usage: land_e0.py --scores DIR --gt GT --out OUT [--ids 200:230 --nboot 30 --seed 0 --ppnorm 1]"""
import os, argparse, json, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--gt', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--ids', default='200:230'); ap.add_argument("--nboot", type=int, default=20); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--ppnorm', type=int, default=1); ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
i0, i1 = [int(x) for x in a.ids.split(':')]; K = i1 - i0 + 1; tg = np.arange(K, dtype=float)
conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int)
S = np.load(a.scores + '/scores.npy'); M = np.load(a.scores + '/meta.npz'); ppn = M['ppnoise']
rows = np.where(np.isfinite(S).all(1) & (conf >= i0) & (conf <= i1))[0]; tt = conf[rows] - i0; N0 = len(rows)
ell = S[rows].astype(np.float64)
if a.ppnorm: ell = ell / np.maximum(ppn[rows], 1e-6)[:, None]
ell -= ell.max(1, keepdims=True)
log('%d swing particles scored against %d truth maps (ppnorm %d)' % (N0, K, a.ppnorm))
def lse(x): m = x.max(1, keepdims=True); return (m + np.log(np.exp(x - m).sum(1, keepdims=True)))[:, 0]
# ---- likelihood calibration (oracle diagnostic): temperature maximising the true-label posterior under the true prior
Nt = np.bincount(tt, minlength=K); lp0 = np.log(Nt / Nt.sum())
def cond_ll(tau): z = ell / tau + lp0[None]; return (z[np.arange(N0), tt] - lse(z)).mean()
taus = np.exp(np.linspace(np.log(0.05), np.log(20.), 81)); cl = np.array([cond_ll(t) for t in taus]); tau_o = float(taus[np.argmax(cl)])
am = ell.argmax(1); post = np.exp(ell - lse(ell)[:, None]); pm = post @ tg
d_end = None
if Nt[0] and Nt[-1]:
    dd = ell[:, 0] - ell[:, -1]; x0, x1 = dd[tt == 0], dd[tt == K - 1]; d_end = float((x0.mean() - x1.mean()) / np.sqrt(0.5 * (x0.var() + x1.var())))
log('calibration: tau* %.3f (cond. log-lik %.4f at tau 1, %.4f at tau*; chance %.4f); argmax exact %.3f, |err| mean %.2f; '
    'flat-prior posterior mean corr %.3f; end-to-end d\' per image %s' % (tau_o, cond_ll(1.), cl.max(), (lp0[tt]).mean(), (am == tt).mean(),
    np.abs(am - tt).mean(), np.corrcoef(pm, tt)[0, 1], 'n/a' if d_end is None else '%.2f' % d_end))
# ---- known landscapes
def target(name):
    if name == 'flat': g = np.zeros(K)
    elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
    elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
    elif name == 'shoulder':
        p = 0.95 * np.exp(-0.5 * ((tg - 20.) / 4.) ** 2) / np.exp(-0.5 * ((tg - 20.) / 4.) ** 2).sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
    else: raise SystemExit('unknown shape ' + name)
    p = np.exp(-g); return p / p.sum()
def subsample(name, seed):
    if name == 'sim': return np.arange(N0)
    p = target(name); s = np.min(Nt / np.maximum(p, 1e-12)); n = np.floor(s * p).astype(int); r = np.random.default_rng(seed)
    return np.concatenate([r.choice(np.where(tt == k)[0], n[k], replace=False) for k in range(K)])
def F(L, p): d = L @ p; return p * (L.T @ (1. / d)) / L.shape[0]
def obj(L, p): return np.log(L @ p).mean()
def npmle(L, p0=None, nit=4000, tol=1e-10):   # EM for the mixing weights on the grid, SQUAREM-accelerated (monotone fallback)
    p = np.full(L.shape[1], 1. / L.shape[1]) if p0 is None else p0.copy(); f = obj(L, p)
    for it in range(nit):
        p1 = F(L, p); p2 = F(L, p1); r = p1 - p; v = p2 - 2. * p1 + p; nv = np.linalg.norm(v)
        if nv < 1e-14: p = p2; break
        al = min(-np.linalg.norm(r) / nv, -1.); pq = np.maximum(p - 2. * al * r + al * al * v, 0.); pq /= pq.sum(); pq = F(L, pq); fq = obj(L, pq)
        if not np.isfinite(fq) or fq < obj(L, p2): pq, fq = p2, obj(L, p2)
        done = np.abs(pq - p).max() < tol; p, f = pq, fq
        if done: break
    return p
def sm3(p): q = np.convolve(p, np.ones(3), 'same') / np.convolve(np.ones_like(p), np.ones(3), 'same'); return q
def metrics(p, pt):
    G = -np.log(np.maximum(sm3(p), 1e-12)); return dict(tv=0.5 * float(np.abs(p - pt).sum()), low=float(p[:15].sum() + 0.5 * p[15]), minor=float(p[:4].sum()),
                                                     barrier=float(G[15] - min(G[5], G[25])))
res = {}; curves = {}
for name in a.shapes.split(','):
    idx = subsample(name, a.seed); pt = np.bincount(tt[idx], minlength=K) / len(idx)
    out = {}
    for tau_name, tau in (('tau1', 1.), ('tauopt', tau_o)):
        L = np.exp(ell[idx] / tau); L /= L.max(1, keepdims=True)
        pe = npmle(L); hard = np.bincount(L.argmax(1), minlength=K) / len(idx)
        pmc = np.clip(np.rint((L / L.sum(1, keepdims=True)) @ tg), 0, K - 1).astype(int); point = np.bincount(pmc, minlength=K) / len(idx)
        rb = np.random.default_rng(a.seed + 7); boots = []
        for b in range(a.nboot):
            ib = rb.integers(0, len(idx), len(idx)); boots.append(npmle(L[ib], p0=pe, nit=400, tol=1e-8))
        boots = np.array(boots); out[tau_name] = dict(npmle=metrics(pe, pt), hard=metrics(hard, pt), point=metrics(point, pt),
                                                      npmle_boot={k: [float(np.percentile([metrics(bb, pt)[k] for bb in boots], q)) for q in (10, 50, 90)] for k in ('tv', 'low', 'minor', 'barrier')})
        curves['%s_%s' % (name, tau_name)] = dict(true=pt, npmle=pe, hard=hard, point=point, lo=np.percentile(boots, 10, 0), hi=np.percentile(boots, 90, 0))
    out['truth'] = metrics(pt, pt); out['n'] = int(len(idx)); res[name] = out
    t1 = out['tau1']; log('%-9s n %6d | true low %.3f minor %.3f barrier %+.2f | NPMLE tv %.3f low %.3f minor %.3f barrier %+.2f [%+.2f %+.2f] | '
                          'hard tv %.3f barrier %+.2f minor %.3f | point tv %.3f barrier %+.2f minor %.3f | tau* NPMLE tv %.3f barrier %+.2f' % (
        name, len(idx), out['truth']['low'], out['truth']['minor'], out['truth']['barrier'], t1['npmle']['tv'], t1['npmle']['low'], t1['npmle']['minor'],
        t1['npmle']['barrier'], t1['npmle_boot']['barrier'][0], t1['npmle_boot']['barrier'][2], t1['hard']['tv'], t1['hard']['barrier'], t1['hard']['minor'],
        t1['point']['tv'], t1['point']['barrier'], t1['point']['minor'], out['tauopt']['npmle']['tv'], out['tauopt']['npmle']['barrier']))
json.dump(dict(tau_opt=tau_o, n_particles=N0, results=res), open(a.out + '/e0.json', 'w'), indent=1)
np.savez(a.out + '/e0_curves.npz', **{k + '__' + kk: vv for k, v in curves.items() for kk, vv in v.items()})
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    names = a.shapes.split(','); fig, ax = plt.subplots(2, len(names), figsize=(3.2 * len(names), 6), squeeze=False)
    for j, name in enumerate(names):
        c = curves[name + '_tau1']
        for i, f in enumerate((lambda p: p, lambda p: -np.log(np.maximum(sm3(p), 1e-6)) + np.log(sm3(c['true']).max()))):
            x = ax[i, j]; x.plot(tg, f(c['true']), 'k-', lw=2, label='truth'); x.plot(tg, f(c['npmle']), 'b-', label='NPMLE')
            x.fill_between(tg, f(c['lo']), f(c['hi']), color='b', alpha=0.2); x.plot(tg, f(c['hard']), color='C1', label='argmax'); x.plot(tg, f(c['point']), color='C2', label='point coord')
            x.set_xlabel('swing t (A of tip travel)'); x.set_ylabel('occupancy' if i == 0 else 'G (kT)'); x.set_title(name)
            if i == 1: x.set_ylim(-0.5, 6)
        if j == 0: ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + '/e0_landscapes.png', dpi=110)
except Exception as e: log('plot skipped: %s' % e)
log('wrote %s/e0.json, e0_curves.npz, e0_landscapes.png' % a.out)
