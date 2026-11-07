#!/usr/bin/env python3
"""land_reweight.py: ensemble reweighting of a path (occupancy p(s) over the K path maps) from land_pscore.py likelihoods, validated on
imposed landscapes. Synthetic mode (--gt): for each score variant (l0 = tree pose, lprof = pose optimised, lmarg = pose marginalised) and
likelihood temperature (1, pixel-split, oracle tau* from the true labels), every imposed shape (subsampling on the true t) is fitted by
the grid NPMLE and by the smooth fit (land_deconv.fit_pen) on a penalty grid; the penalty is CHOSEN ON TRUTH (lowest mean TV over the
shapes, barrier error as tie-break) and bootstrap bands are given at that choice. Real mode (--lam, --variant, --tau): p(s) with bands,
G(s), each particle's posterior mean s, and window labels (master rows) for reconstructions along the path.
usage: land_reweight.py --scores NPZ --out DIR [--gt GT --ids 200:230 --shapes ...] | [--real --variant lmarg --tau split --lam 1e-3 --nm NM]"""
import os, sys, json, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--out', required=True); ap.add_argument('--gt', default=''); ap.add_argument('--ids', default='200:230')
ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder'); ap.add_argument('--variants', default='l0,lprof,lmarg'); ap.add_argument('--taus', default='1,split,opt')
ap.add_argument('--lams', default='0,1e-6,1e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1'); ap.add_argument('--nboot', type=int, default=20); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--real', action='store_true'); ap.add_argument('--variant', default='lmarg'); ap.add_argument('--tau', default='split'); ap.add_argument('--lam', type=float, default=1e-3)
ap.add_argument('--nm', type=int, default=0); ap.add_argument('--nwin', type=int, default=8); ap.add_argument('--exclude', default='', help='real mode: master-row flags (nonzero = drop)')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
Z = np.load(a.scores); rows = Z['rows']; K = Z['l0'].shape[1]; tg = np.arange(K, dtype=float); tau_split = float(np.median(Z['tau_split'])); half = Z['half'].astype(int); eo = Z['eo'].astype(int)
cvh = eo if ((eo == 0).any() and (eo == 1).any()) else half if (half.max() > 0) else np.random.default_rng(a.seed).integers(0, 2, len(rows))
log('%s: %d particles x %d path maps; pixel-split temperature %.3f (pairs %s)' % (a.scores, len(rows), K, tau_split, ' '.join('%.3f' % x for x in Z['tau_split'])))
def lik(var, tau, idx=slice(None)):
    l = Z[var][idx].astype(np.float64) / tau; l -= l.max(1, keepdims=True); return np.exp(l)
def lse(x): m = x.max(1, keepdims=True); return (m + np.log(np.exp(x - m).sum(1, keepdims=True)))[:, 0]
lams = [float(x) for x in a.lams.split(',')]
if not a.real:
    i0, i1 = [int(x) for x in a.ids.split(':')]; conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int); tt = conf[rows] - i0; ok = (tt >= 0) & (tt < K); tt = np.where(ok, tt, -1)
    Nt = np.bincount(tt[ok], minlength=K); lp0 = np.log(np.maximum(Nt, 1) / Nt.sum())
    def tau_opt(var):
        l = Z[var][ok].astype(np.float64); taus = np.exp(np.linspace(np.log(0.1), np.log(20.), 61)); best = None
        for t in taus:
            zz = l / t + lp0[None]; v = (zz[np.arange(ok.sum()), tt[ok]] - lse(zz)).mean()
            if best is None or v > best[1]: best = (t, v)
        return float(best[0])
    def target(name):
        if name == 'flat': g = np.zeros(K)
        elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
        elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
        elif name == 'shoulder':
            m_ = np.exp(-0.5 * ((tg - 20.) / 4.) ** 2); p = 0.95 * m_ / m_.sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
        p = np.exp(-g); return p / p.sum()
    def subsample(name):
        cod = np.where(ok)[0]
        if name == 'sim': return cod
        p = target(name); s = np.min(Nt / np.maximum(p, 1e-12)); nk = np.floor(s * p).astype(int); r = np.random.default_rng(a.seed)
        return np.concatenate([r.choice(cod[tt[cod] == k], nk[k], replace=False) for k in range(K)])
    shapes = a.shapes.split(','); SUB = {s: subsample(s) for s in shapes}; PT = {s: np.bincount(tt[SUB[s]], minlength=K) / len(SUB[s]) for s in shapes}
    res = dict(tau_split=tau_split, table={}, chosen={}); curves = {}
    for var in a.variants.split(','):
        to = tau_opt(var); log('variant %s: oracle tau* %.3f, pixel-split %.3f; argmax exact %.3f' % (var, to, tau_split, (Z[var][ok].argmax(1) == tt[ok]).mean()))
        for tn in a.taus.split(','):
            tau = 1. if tn == '1' else tau_split if tn == 'split' else to; key = '%s_tau%s' % (var, tn); tab = {}
            for s in shapes:
                L = lik(var, tau, SUB[s]); pt = PT[s]; row = {}
                for lm in lams:
                    p = LD.npmle(L) if lm == 0 else LD.fit_pen(L, lm)[1]; row[str(lm)] = dict(LD.metrics(p, pt), p=p.tolist())
                tab[s] = row
            mt = {lm: np.mean([tab[s][str(lm)]['tv'] for s in shapes]) for lm in lams}
            mb = {lm: np.mean([abs(tab[s][str(lm)]['barrier'] - LD.metrics(PT[s], PT[s])['barrier']) for s in shapes if s.startswith('dw') or s in ('flat', 'broad')]) for lm in lams}
            lb = min(lams, key=lambda x: (round(mt[x], 3), mb[x])); res['table'][key] = {s: {lm: {k_: v_ for k_, v_ in tab[s][lm].items() if k_ != 'p'} for lm in tab[s]} for s in shapes}
            res['chosen'][key] = dict(lam=lb, mean_tv=mt[lb], mean_barrier_err=mb[lb], tau=tau)
            log('  %-14s tau %.3f: mean TV per lam %s | chosen lam %g (mean TV %.3f, mean |barrier err| %.2f kT)' % (key, tau, ' '.join('%g:%.3f' % (x, mt[x]) for x in lams), lb, mt[lb], mb[lb]))
            for s in shapes:
                L = lik(var, tau, SUB[s]); p = LD.npmle(L) if lb == 0 else LD.fit_pen(L, lb)[1]; B = []; rng = np.random.default_rng(a.seed + 3)
                for b in range(a.nboot):
                    ib = rng.integers(0, len(L), len(L)); B.append(LD.npmle(L[ib], p0=p, nit=400, tol=1e-8) if lb == 0 else LD.fit_pen(L[ib], lb, G0=-np.log(np.maximum(p, 1e-12)))[1])
                m_, bd = LD.metrics(p, PT[s]), LD.band(np.array(B), PT[s]); t_ = LD.metrics(PT[s], PT[s]); curves['%s__%s' % (key, s)] = np.c_[PT[s], p, np.percentile(B, 10, 0), np.percentile(B, 90, 0)]
                res['chosen'][key][s] = dict(m=m_, band=bd, truth=t_)
                log('    %-9s TRUE low %.3f minor %.3f barrier %+.2f | est tv %.3f low %.3f minor %.3f [%.3f %.3f] barrier %+.2f [%+.2f %+.2f]' % (s, t_['low'], t_['minor'], t_['barrier'],
                    m_['tv'], m_['low'], m_['minor'], bd['minor'][0], bd['minor'][2], m_['barrier'], bd['barrier'][0], bd['barrier'][2]))
    json.dump(res, open(a.out + '/reweight.json', 'w'), indent=1); np.savez(a.out + '/reweight_curves.npz', **curves)
    best = min(res['chosen'], key=lambda k_: (res['chosen'][k_]['mean_tv'], res['chosen'][k_]['mean_barrier_err']))
    log('RESULT best setting %s (lam %g): mean TV %.3f, mean |barrier err| %.2f kT' % (best, res['chosen'][best]['lam'], res['chosen'][best]['mean_tv'], res['chosen'][best]['mean_barrier_err']))
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        keys = [k_ for k_ in res['chosen']]; fig, ax = plt.subplots(len(keys), len(shapes), figsize=(2.9 * len(shapes), 2.3 * len(keys)), squeeze=False)
        for i, k_ in enumerate(keys):
            for j, s in enumerate(shapes):
                cv = curves['%s__%s' % (k_, s)]; ref = np.log(LD.sm3(cv[:, 0]).max()); x = ax[i, j]; G_ = lambda p: -np.log(np.maximum(LD.sm3(p), 1e-6)) + ref
                x.plot(tg, G_(cv[:, 0]), 'k-', lw=2); x.plot(tg, G_(cv[:, 1]), 'b-'); x.fill_between(tg, G_(cv[:, 2]), G_(cv[:, 3]), color='b', alpha=0.2); x.set_ylim(-0.5, 6)
                if j == 0: x.set_ylabel('%s\nG (kT)' % k_, fontsize=7)
                if i == 0: x.set_title(s)
        fig.tight_layout(); fig.savefig(a.out + '/reweight_G.png', dpi=100)
    except Exception as e: log('plot skipped: %s' % e)
else:
    sel = np.arange(len(rows))
    if a.exclude: sel = np.where(np.loadtxt(a.exclude).astype(int)[rows] == 0)[0]; log('excluded %d of %d particles (%s)' % (len(rows) - len(sel), len(rows), a.exclude))
    rows, cvh = rows[sel], cvh[sel]
    tau = 1. if a.tau == '1' else tau_split if a.tau == 'split' else float(a.tau); L = lik(a.variant, tau, sel); p = LD.npmle(L) if a.lam == 0 else LD.fit_pen(L, a.lam)[1]
    B = []; rng = np.random.default_rng(a.seed)
    for b in range(a.nboot):
        ib = rng.integers(0, len(L), len(L)); B.append(LD.npmle(L[ib], p0=p, nit=400, tol=1e-8) if a.lam == 0 else LD.fit_pen(L[ib], a.lam, G0=-np.log(np.maximum(p, 1e-12)))[1])
    B = np.array(B); post = L * p[None]; post /= post.sum(1, keepdims=True); sm = post @ (tg / (K - 1.))
    G_ = -np.log(np.maximum(p, 1e-12)); G_ -= G_.min(); Gb = -np.log(np.maximum(B, 1e-12)); Gb -= Gb.min(1, keepdims=True)
    # half-stack replicate: the fit on each half separately
    ph = [LD.npmle(L[cvh == h]) if a.lam == 0 else LD.fit_pen(L[cvh == h], a.lam)[1] for h in (0, 1)]
    log('variant %s, tau %.3f, lam %g: p(s) over %d path maps (s = 0 first map .. 1 last map)' % (a.variant, tau, a.lam, K))
    for k in range(K): log('  s %.3f  p %.4f [%.4f %.4f]  halves %.4f / %.4f  G %.2f [%.2f %.2f] kT' % (k / (K - 1.), p[k], *np.percentile(B[:, k], [10, 90]), ph[0][k], ph[1][k], G_[k], *np.percentile(Gb[:, k], [10, 90])))
    qs = np.percentile(sm, np.linspace(0, 100, a.nwin + 1)); wl = np.clip(np.searchsorted(qs, sm, side='right') - 1, 0, a.nwin - 1) + 1
    if a.nm: lab = np.zeros(a.nm, int); lab[rows] = wl; np.savetxt(a.out + '/windows_path.txt', lab, fmt='%d')
    json.dump(dict(variant=a.variant, tau=tau, lam=a.lam, p=p.tolist(), lo=np.percentile(B, 10, 0).tolist(), hi=np.percentile(B, 90, 0).tolist(), halves=[x.tolist() for x in ph],
                   G=G_.tolist(), s_mean_quantiles=qs.tolist()), open(a.out + '/reweight_real.json', 'w'), indent=1)
    np.savez(a.out + '/posterior.npz', rows=rows, smean=sm, post=post.astype(np.float32)); log('RESULT real: wrote %s; DONE' % a.out)
