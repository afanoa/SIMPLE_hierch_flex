#!/usr/bin/env python3
"""Particle-level comparator from a subcls_embed run's npz: the bin-centred node latent projected on the most split-half-reliable
direction (truth-free, chosen with the sub-classes), binned PER PARTICLE into K children by (width) equal-width bins and (binq)
view-balanced within-bin quantiles; per rule the swing eta2 (truth) and the children's view balance (KL of the bin histogram vs the
node's, and for a random subset of the same size). Writes children_pp_<rule>.txt per row.
usage: subcls_ppbins.py EMBED_DIR NROWS GT [--codes swing:200:230] [--nchild 10] [--dir 1]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('nrows', type=int); ap.add_argument('gt'); ap.add_argument('--codes', default='swing:200:230')
ap.add_argument('--nchild', type=int, default=10); ap.add_argument('--dir', type=int, default=1); a = ap.parse_args(); rng = np.random.default_rng(0)
E = np.load(os.path.join(a.d, 'subcls_embed.npz')); p, sid, sbin, Z, V = E['p'], E['sid'], E['sbin'], E['Z'], E['bincentred_V']
binn = sbin[sid]; ub, bi = np.unique(binn, return_inverse=True); nb = len(ub); K = a.nchild
bm = np.stack([np.bincount(bi, Z[:, j], minlength=nb) for j in range(Z.shape[1])], 1) / np.bincount(bi, minlength=nb)[:, None]; Zc = Z - bm[bi]
y = (Zc - Zc.mean(0)) @ V[:, a.dir - 1]
G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(a.nrows + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
nm, lo, hi = a.codes.split(':'); v = np.where((conf >= int(lo)) & (conf <= int(hi)), conf - int(lo), np.nan); dv = np.isfinite(v)
h0 = np.bincount(bi, minlength=nb) / len(bi)
def kl(m):
    h = np.bincount(bi[m], minlength=nb) / m.sum(); r = rng.choice(len(bi), m.sum(), replace=False); hr = np.bincount(bi[r], minlength=nb) / m.sum()
    f = lambda q: (q[q > 0] * np.log(q[q > 0] / h0[q > 0])).sum(); return f(h), f(hr)
def eta2(lab):
    u, inv = np.unique(lab[dv], return_inverse=True); m = np.bincount(inv, v[dv]) / np.bincount(inv); return 1 - ((v[dv] - m[inv]) ** 2).sum() / ((v[dv] - v[dv].mean()) ** 2).sum()
l1, h1 = np.percentile(y, [1, 99]); rules = {'width': np.clip(np.digitize(y, np.linspace(l1, h1, K + 1)[1:-1]), 0, K - 1) + 1}
q = np.zeros(len(y), int)
for b in range(nb):
    m = np.where(bi == b)[0]; r = np.argsort(np.argsort(y[m])); q[m] = np.minimum(r * K // len(m), K - 1) + 1
rules['binq'] = q
print('direction %d of the bin-centred latent, per particle: |corr| with %s %.3f' % (a.dir, nm, abs(np.corrcoef(y[dv], v[dv])[0, 1])))
for rn, lab in rules.items():
    out = np.zeros(a.nrows, int); out[p - 1] = lab; np.savetxt(os.path.join(a.d, 'children_pp_%s.txt' % rn), out, fmt='%d')
    ks = [kl(lab == c) for c in range(1, K + 1)]
    print('[per particle, %s] %s eta2 %.3f; view KL per child %s (random subsets %.3f-%.3f); per child mean %s: %s'
          % (rn, nm, eta2(lab), ' '.join('%.3f' % k_[0] for k_ in ks), min(k_[1] for k_ in ks), max(k_[1] for k_ in ks), nm,
             ' '.join('%.1f' % np.nanmean(v[dv & (lab == c)]) for c in range(1, K + 1))))
