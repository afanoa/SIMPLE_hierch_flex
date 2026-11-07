#!/usr/bin/env python3
"""Post-analysis of a subcls_embed2.py run (reads its subcls_embed.npz; nothing is recomputed upstream): children of the leading
reproducible coordinate under two rules, (width) equal-width bins over sub-classes (the embed2 default) and (binq) VIEW-BALANCED
within-bin quantiles (inside every direction bin the sub-classes are ordered along the coordinate and split by the member-weighted
CDF, so every child takes ~1/K of every bin), with per child: members, view balance (KL of its bin histogram vs the node's, and
the same for a random subset of equal size), and truth (--codes as in embed2). Writes children_<tag>_<rule>.txt (per row).
usage: subcls_post.py EMBED_DIR NROWS [--tag bincentred] [--nchild 10] [--dir 1] [--gt gt --codes swing:200:230]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('nrows', type=int); ap.add_argument('--tag', default='bincentred')
ap.add_argument('--nchild', type=int, default=10); ap.add_argument('--dir', type=int, default=1); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='')
a = ap.parse_args(); E = np.load(os.path.join(a.d, 'subcls_embed.npz')); rng = np.random.default_rng(0)
p, sid, nS, sbin = E['p'], E['sid'], E['nS'], E['sbin']; y = E[a.tag + '_Ys'][:, a.dir - 1]; S = len(nS); K = a.nchild; binn = sbin[sid]
ub, bi = np.unique(binn, return_inverse=True); nb = len(ub); h0 = np.bincount(bi, minlength=nb) / len(bi)
def kl(lab, c):
    m = lab == c; h = np.bincount(bi[m], minlength=nb) / max(m.sum(), 1); r = rng.choice(len(bi), m.sum(), replace=False); hr = np.bincount(bi[r], minlength=nb) / max(m.sum(), 1)
    f = lambda q: (q[q > 0] * np.log(q[q > 0] / h0[q > 0])).sum(); return f(h), f(hr)
rules = {}
o = np.argsort(y); cw = np.cumsum(nS[o]) / nS.sum(); lo, hi = y[o][np.searchsorted(cw, 0.01)], y[o][np.searchsorted(cw, 0.99)]
rules['width'] = (np.clip(np.digitize(y, np.linspace(lo, hi, K + 1)[1:-1]), 0, K - 1) + 1)[sid]
cs = np.zeros(S, int)
for b in np.unique(sbin):
    m = np.where(sbin == b)[0]; m = m[np.argsort(y[m])]; c = (np.cumsum(nS[m]) - 0.5 * nS[m]) / nS[m].sum(); cs[m] = np.minimum((c * K).astype(int), K - 1) + 1
rules['binq'] = cs[sid]
TV = {}
if a.gt:
    G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(a.nrows + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
    for c in [c for c in a.codes.split(',') if c]:
        nm, l0, h0_ = c.split(':'); TV[nm] = np.where((conf >= int(l0)) & (conf <= int(h0_)), conf - int(l0), np.nan)
def eta2(v, lab):
    u, inv = np.unique(lab, return_inverse=True); m = np.bincount(inv, v) / np.bincount(inv); return 1 - ((v - m[inv]) ** 2).sum() / ((v - v.mean()) ** 2).sum()
for rn, lab in rules.items():
    out = np.zeros(a.nrows, int); out[p - 1] = lab; np.savetxt(os.path.join(a.d, 'children_%s_%s.txt' % (a.tag, rn)), out, fmt='%d')
    kls = [kl(lab, c) for c in range(1, K + 1)]
    line = '[%s %s] children (members, view KL vs random-subset KL):' % (a.tag, rn) + ''.join(' %d:%d(%.3f/%.3f)' % (c, (lab == c).sum(), *kls[c - 1]) for c in range(1, K + 1))
    print(line)
    for nm, v in TV.items():
        dv = np.isfinite(v)
        print('   TRUTH %s eta2 %.3f; per child mean (share): %s' % (nm, eta2(v[dv], lab[dv]), ' '.join('%.1f(%.0f%%)' % (np.nanmean(v[dv & (lab == c)]) if (dv & (lab == c)).any() else np.nan, 100 * dv[lab == c].mean()) for c in range(1, K + 1))))
