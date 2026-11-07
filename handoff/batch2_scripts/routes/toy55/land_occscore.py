#!/usr/bin/env python3
"""land_occscore.py: scores a synthetic land_occ.py run in the READOUT units the estimate lives in. The window maps' closed-site readout is a
compressed linear image of the particles' true occupancy (reconstruction and pose blur), so the truth is mapped through the windows' fit
readout = alpha + beta * true mean occupancy (truth used for scoring only). Reports, against the mapped truth: TV, the masses of the two
halves of the occupancy range of the truth, the minor (highest-occupancy 5 % of the truth) mass, the dip depth between the truth's two
largest modes (kT, 3-bin smoothed), and the oracle noise scale in these units (residual of lam about E[lam | true occupancy] x Fisher).
usage: land_occscore.py RUNDIR TRUTHOCC (RUNDIR holds e1b_particles.npz and occ/occ.json; TRUTHOCC = land_truthocc.py output)"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
run = sys.argv[1]; R = json.load(open(run + '/occ/occ.json')); Z = np.load(run + '/e1b_particles.npz')
og = np.array(R['grid']); W = R['windows']; tw = np.array([w['true'] for w in W]); rw = np.array([np.mean(w['occ']) for w in W])
beta, alpha = np.polyfit(tw, rw, 1); print('window readout = %.3f + %.3f x true occupancy (r %.3f)' % (alpha, beta, np.corrcoef(tw, rw)[0, 1]))
occ_tab = np.array(open(sys.argv[2]).read().split(), float)
t = Z['t']; occ = np.where(t >= 0, occ_tab[np.maximum(t, 0)], np.nan); m = np.isfinite(occ); mapped = alpha + beta * occ[m]
step = og[1] - og[0]; pt = np.bincount(np.clip(np.rint((mapped - og[0]) / step).astype(int), 0, len(og) - 1), minlength=len(og)) / m.sum()
def sm3(p): return LD.sm3(np.asarray(p))
lo, hi = np.percentile(mapped, [5, 95]); mid = 0.5 * (lo + hi); top5 = np.percentile(mapped, 95)
def dip(p):   # depth of the deepest minimum between the two largest maxima of the 3-bin smoothed curve (kT); 0 if unimodal
    q = sm3(p); G = -np.log(np.maximum(q, 1e-12)); mx = [i for i in range(1, len(q) - 1) if q[i] >= q[i - 1] and q[i] >= q[i + 1] and q[i] > 1e-4]
    if len(mx) < 2: return 0.
    mx = sorted(sorted(mx, key=lambda i: -q[i])[:2]); a_, b_ = mx; j = a_ + int(np.argmax(G[a_:b_ + 1])); return float(G[j] - max(G[a_], G[b_]))
def score(p): p = np.asarray(p); return dict(tv=0.5 * float(np.abs(p - pt).sum()), upper=float(p[og >= mid].sum()), minor=float(p[og >= top5].sum()), dip=dip(p))
lam = Z['lam']; I = Z['I']; lamf = (Z['lamA'] * Z['IA'] + Z['lamB'] * Z['IB']) / np.maximum(Z['IA'] + Z['IB'], 1e-30); If = Z['IA'] + Z['IB']
sg = R.get('orient', 1.); lamf = sg * lamf; ks = np.unique(t[t >= 0]); El = {k: lamf[t == k].mean() for k in ks}; res = lamf[t >= 0] - np.array([El[k] for k in t[t >= 0]])
c_or = float(np.mean(res ** 2 * If[t >= 0])); print('noise scale: c_reg %.3f, oracle %.3f (residual about E[lam | t] x Fisher)' % (R['c_reg'], c_or))
out = dict(alpha=alpha, beta=beta, truth=score(pt), c_oracle=c_or)
for vn in ('c_reg', 'c_split'):
    if vn in R: out[vn] = score(R[vn]['p']); print('%-8s %s' % (vn, ' '.join('%s %.3f' % kv for kv in out[vn].items())))
print('truth    %s' % ' '.join('%s %.3f' % kv for kv in out['truth'].items()))
json.dump(out, open(run + '/occ/occscore.json', 'w'), indent=1)
