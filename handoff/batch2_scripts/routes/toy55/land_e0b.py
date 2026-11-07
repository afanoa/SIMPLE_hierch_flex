#!/usr/bin/env python3
"""land_e0b.py: re-analysis of the oracle scores (land_e0.sbatch: swing particles vs the true swing maps) with the SMOOTH cross-half
deconvolution (land_deconv.fit_cv_boot) and a truth-free likelihood temperature from the pixel split: tau_split = variance of the
set-A minus set-B score differences of map pairs over the model's prediction (the score noise covariances), median over pairs.
Estimators per imposed landscape: argmax histogram, point coordinate (flat-prior posterior mean), NPMLE, smooth CV fit at tau 1 / tau_split /
the oracle tau*. usage: land_e0b.py --scores DIR --gt GT --out OUT [--shapes ... --nboot 20]"""
import os, sys, argparse, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--gt', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--ids', default='200:230'); ap.add_argument('--nboot', type=int, default=20); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--shapes', default='sim,flat,broad,dw2,dw4,shoulder'); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
i0, i1 = [int(x) for x in a.ids.split(':')]; K = i1 - i0 + 1; tg = np.arange(K, dtype=float)
conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int)
SA = np.load(a.scores + '/scores_pA.npy'); SB = np.load(a.scores + '/scores_pB.npy'); M = np.load(a.scores + '/meta.npz'); ppn = M['ppnoise']; EO = M['eo']
rows = np.where(np.isfinite(SA).all(1) & np.isfinite(SB).all(1) & (conf >= i0) & (conf <= i1))[0]; tt = conf[rows] - i0; N0 = len(rows)
pn = np.maximum(ppn[rows], 1e-6)[:, None]; A_ = SA[rows] / pn; B_ = SB[rows] / pn; ell = A_ + B_; ell -= ell.max(1, keepdims=True); eo = EO[rows]
half = eo.astype(int) if (eo == 0).any() and (eo == 1).any() else np.random.default_rng(a.seed).integers(0, 2, N0)
# truth-free temperature: split differences of map-pair score contrasts vs the model's noise prediction
GA = {e: M['gramA%d' % e] for e in (0, 1)}; GB = {e: M['gramB%d' % e] for e in (0, 1)}; rat = []
gsc = np.mean(1. / pn[:, 0])
for j in range(0, K, 3):
    for k in range(j + 3, K, 3):
        dA = A_[:, j] - A_[:, k]; dB = B_[:, j] - B_[:, k]; emp = np.var(dA - dB)
        pred = np.mean([(GA[e][j, j] + GA[e][k, k] - 2 * GA[e][j, k] + GB[e][j, j] + GB[e][k, k] - 2 * GB[e][j, k]) for e in (0, 1) if np.trace(GA[e]) > 0]) * gsc
        rat.append(emp / max(pred, 1e-30))
tau_split = float(np.median(rat)); log('pixel-split temperature: median %.3f (pairs %d, 10/90 %% %.3f/%.3f)' % (tau_split, len(rat), *np.percentile(rat, [10, 90])))
Nt = np.bincount(tt, minlength=K); lp0 = np.log(Nt / Nt.sum())
def lse(x): m = x.max(1, keepdims=True); return (m + np.log(np.exp(x - m).sum(1, keepdims=True)))[:, 0]
def cond_ll(tau): z = ell / tau + lp0[None]; return (z[np.arange(N0), tt] - lse(z)).mean()
taus = np.exp(np.linspace(np.log(0.05), np.log(20.), 81)); tau_o = float(taus[np.argmax([cond_ll(t) for t in taus])])
log('oracle temperature tau* %.3f; %d particles, halves %d / %d' % (tau_o, N0, (half == 0).sum(), (half == 1).sum()))
# gradient check of the penalised objective (small random problem)
_L = np.random.default_rng(1).random((200, 7)); _G = np.random.default_rng(2).normal(size=7); _lam = 0.01
D2 = np.diff(np.eye(7), 2, axis=0); _P = D2.T @ D2
def _f(G): p = LD._pk(G); return -np.log(_L @ p).mean() + _lam * G @ _P @ G
num = np.array([(_f(_G + 1e-6 * np.eye(7)[k]) - _f(_G - 1e-6 * np.eye(7)[k])) / 2e-6 for k in range(7)])
p_ = LD._pk(_G); r_ = (_L * p_) / (_L @ p_)[:, None]; ana = (r_.mean(0) - p_) + 2 * _lam * _P @ _G; log('gradient check: max |num - analytic| %.2e' % np.abs(num - ana).max())
def target(name):
    if name == 'flat': g = np.zeros(K)
    elif name.startswith('dw'): h = float(name[2:]); g = h * (((tg - 15.) / 10.) ** 2 - 1.) ** 2
    elif name == 'broad': g = 0.5 * ((tg - 15.) / 7.) ** 2
    elif name == 'shoulder':
        m_ = np.exp(-0.5 * ((tg - 20.) / 4.) ** 2); p = 0.95 * m_ / m_.sum() + 0.05 * (tg <= 3) / (tg <= 3).sum(); return p / p.sum()
    p = np.exp(-g); return p / p.sum()
def subsample(name, seed):
    if name == 'sim': return np.arange(N0)
    p = target(name); s = np.min(Nt / np.maximum(p, 1e-12)); n = np.floor(s * p).astype(int); r = np.random.default_rng(seed)
    return np.concatenate([r.choice(np.where(tt == k)[0], n[k], replace=False) for k in range(K)])
res = dict(tau_split=tau_split, tau_opt=tau_o, n=N0); curves = {}
for name in a.shapes.split(','):
    idx = subsample(name, a.seed); pt = np.bincount(tt[idx], minlength=K) / len(idx); out = dict(n=int(len(idx)), truth=LD.metrics(pt, pt))
    L1 = np.exp(ell[idx]); hard = np.bincount(L1.argmax(1), minlength=K) / len(idx)
    point = np.bincount(np.clip(np.rint((L1 / L1.sum(1, keepdims=True)) @ tg), 0, K - 1).astype(int), minlength=K) / len(idx)
    out['hard'] = LD.metrics(hard, pt); out['point'] = LD.metrics(point, pt); out['npmle'] = LD.metrics(LD.npmle(L1), pt); curves[name] = dict(true=pt, hard=hard, point=point)
    for tn, tau in (('tau1', 1.), ('tausplit', tau_split), ('tauopt', tau_o)):
        L = np.exp(ell[idx] / tau); fr = LD.fit_cv_boot(L, half[idx], nboot=a.nboot, seed=a.seed)
        out['pen_' + tn] = dict(m=LD.metrics(fr['p'], pt), band=LD.band(fr['boot'], pt), lam=fr['lam'], cv=fr['cvtab'].tolist())
        curves[name]['pen_' + tn] = fr['p']; curves[name]['pen_' + tn + '_lo'] = np.percentile(fr['boot'], 10, 0); curves[name]['pen_' + tn + '_hi'] = np.percentile(fr['boot'], 90, 0)
    res[name] = out; t_ = out['truth']
    log('%-9s n %6d | TRUE low %.3f minor %.3f barrier %+.2f' % (name, len(idx), t_['low'], t_['minor'], t_['barrier']))
    for k_ in ('hard', 'point', 'npmle'): m_ = out[k_]; log('    %-13s tv %.3f low %.3f minor %.3f barrier %+.2f' % (k_, m_['tv'], m_['low'], m_['minor'], m_['barrier']))
    for tn in ('tau1', 'tausplit', 'tauopt'):
        m_, b_ = out['pen_' + tn]['m'], out['pen_' + tn]['band']; log('    pen_%-9s tv %.3f low %.3f minor %.3f [%.3f %.3f] barrier %+.2f [%+.2f %+.2f] (lam %.0e)' % (tn, m_['tv'], m_['low'], m_['minor'],
            b_['minor'][0], b_['minor'][2], m_['barrier'], b_['barrier'][0], b_['barrier'][2], out['pen_' + tn]['lam']))
json.dump(res, open(a.out + '/e0b.json', 'w'), indent=1); np.savez(a.out + '/e0b_curves.npz', **{k + '__' + kk: vv for k, v in curves.items() for kk, vv in v.items()})
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    names = a.shapes.split(','); fig, ax = plt.subplots(2, len(names), figsize=(3.3 * len(names), 6.4), squeeze=False)
    for j, name in enumerate(names):
        c = curves[name]; ref = np.log(LD.sm3(c['true']).max())
        for i in (0, 1):
            f = (lambda p: p) if i == 0 else (lambda p: -np.log(np.maximum(LD.sm3(p), 1e-6)) + ref); x = ax[i, j]
            x.plot(tg, f(c['true']), 'k-', lw=2.2, label='truth'); x.plot(tg, f(c['point']), color='C2', label='point coordinate'); x.plot(tg, f(c['hard']), color='C1', label='argmax')
            x.plot(tg, f(c['pen_tausplit']), 'b-', label='deconvolved (split tau)'); x.fill_between(tg, f(c['pen_tausplit_lo']), f(c['pen_tausplit_hi']), color='b', alpha=0.2)
            x.set_xlabel('swing t (A tip travel)'); x.set_ylabel('occupancy' if i == 0 else 'G (kT)'); x.set_title(name)
            if i == 1: x.set_ylim(-0.5, 6)
        if j == 0: ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + '/e0b_landscapes.png', dpi=110)
except Exception as e: log('plot skipped: %s' % e)
log('DONE %s' % a.out)
