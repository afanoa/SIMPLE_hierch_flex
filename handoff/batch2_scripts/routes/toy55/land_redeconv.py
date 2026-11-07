#!/usr/bin/env python3
"""land_redeconv.py: re-run the occupancy deconvolution of a land_occ.py run under several smoothness rules (no maps needed): per-particle
cross-half held-out log-likelihoods for every penalty lam; rules 'max' (best mean held-out log-likelihood: picks lam = 0 on spiky
solutions), '1se' (the largest lam within one standard error of the best), and fixed lam values. Synthetic runs are scored in the readout
units (truth mapped through the windows' readout = alpha + beta * true occupancy, scoring only).
usage: land_redeconv.py RUNDIR [--truthocc TABLE] [--kind e1b|real] [--nboot 20]"""
import os, sys, json, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--truthocc', default=''); ap.add_argument('--kind', default='e1b'); ap.add_argument('--nboot', type=int, default=20)
ap.add_argument('--lams', default='0,1e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.3,1'); a = ap.parse_args()
R = json.load(open(a.run + '/occ/occ.json')); Z = np.load(a.run + ('/e1b_particles.npz' if a.kind == 'e1b' else '/real_particles.npz'))
og = np.array(R['grid']); lam_o = np.array(R['lam_occ']); c_reg = R['c_reg']; sg = R.get('orient', 1.)
half = (Z['eo'] if a.kind == 'e1b' else Z['half']).astype(int); IA, IB = Z['IA'], Z['IB']; If = IA + IB; lamf = sg * (Z['lamA'] * IA + Z['lamB'] * IB) / np.maximum(If, 1e-30)
dl = lamf[:, None] - lam_o[None]; L = np.exp(-0.5 * dl ** 2 * (If / c_reg)[:, None]); L /= np.maximum(L.max(1, keepdims=True), 1e-300)
lams = [float(x) for x in a.lams.split(',')]; HL = np.zeros((len(lams), len(L)))
for i, lm in enumerate(lams):
    for h in (0, 1):
        _, p = LD.fit_pen(L[half == h], lm); HL[i, half != h] = np.log(np.maximum(L[half != h] @ p, 1e-300))
mu = HL.mean(1); se = HL.std(1) / np.sqrt(HL.shape[1]); dse = np.array([np.std(HL[i] - HL[np.argmax(mu)]) / np.sqrt(HL.shape[1]) for i in range(len(lams))])
ib = int(np.argmax(mu)); i1 = max([i for i in range(len(lams)) if mu[i] >= mu[ib] - dse[i]] or [ib])
print('held-out log-lik per particle (minus best): ' + ' '.join('%g:%+.5f(%.5f)' % (lm, m - mu[ib], d) for lm, m, d in zip(lams, mu, dse)))
print('rule max -> lam %g; rule 1se (paired SE) -> lam %g' % (lams[ib], lams[i1]))
pt = None
if a.truthocc:
    W = R['windows']; tw = np.array([w['true'] for w in W]); rw = np.array([np.mean(w['occ']) for w in W]); beta, alpha = np.polyfit(tw, rw, 1)
    tab = np.array(open(a.truthocc).read().split(), float); t = Z['t']; m = t >= 0; mapped = alpha + beta * tab[t[m]]
    pt = np.bincount(np.clip(np.rint((mapped - og[0]) / (og[1] - og[0])).astype(int), 0, len(og) - 1), minlength=len(og)) / m.sum()
    print('truth in readout units (readout = %.3f + %.3f x true): %s' % (alpha, beta, ' '.join('%.3f' % x for x in pt)))
def dip(p):
    q = LD.sm3(np.asarray(p)); G = -np.log(np.maximum(q, 1e-12)); mx = [i for i in range(1, len(q) - 1) if q[i] >= q[i - 1] and q[i] >= q[i + 1] and q[i] > 1e-3]
    if len(mx) < 2: return 0.
    a_, b_ = sorted(sorted(mx, key=lambda i: -q[i])[:2]); j = a_ + int(np.argmax(G[a_:b_ + 1])); return float(G[j] - max(G[a_], G[b_]))
out = dict(lams=lams, heldout=mu.tolist(), dse=dse.tolist(), lam_max=lams[ib], lam_1se=lams[i1], fits={})
rng = np.random.default_rng(0)
for i, lm in enumerate(lams):
    G, p = LD.fit_pen(L, lm); B = []
    if i in (ib, i1) or lm in (1e-3, 1e-2):
        for b in range(a.nboot): ix = rng.integers(0, len(L), len(L)); B.append(LD.fit_pen(L[ix], lm, G0=G)[1])
    row = dict(p=p.tolist(), dip=dip(p))
    if B: row['dip_band'] = [float(np.percentile([dip(b) for b in B], q)) for q in (10, 50, 90)]; row['lo'] = np.percentile(B, 10, 0).tolist(); row['hi'] = np.percentile(B, 90, 0).tolist()
    if pt is not None: row['tv'] = 0.5 * float(np.abs(p - pt).sum())
    out['fits'][str(lm)] = row
    print('lam %-6g dip %.2f%s%s | p %s' % (lm, row['dip'], '' if not B else ' [%.2f %.2f]' % (row['dip_band'][0], row['dip_band'][2]), '' if pt is None else ' tv %.3f' % row['tv'],
          ' '.join('%.3f' % x for x in p)))
if pt is not None: out['truth'] = pt.tolist(); out['truth_dip'] = dip(pt); print('truth dip %.2f' % dip(pt))
json.dump(out, open(a.run + '/occ/redeconv.json', 'w'), indent=1)
