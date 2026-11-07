#!/usr/bin/env python3
"""land_fig_real.py: occupancy p(s) and G(s) along named paths for real reweighting runs (reweight_real.json), with bootstrap bands and the
two half-stack fits. usage: land_fig_real.py OUT.png LABEL:DIR [LABEL:DIR ...]"""
import sys, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out = sys.argv[1]; runs = [x.split(':', 1) for x in sys.argv[2:]]; fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for j, (lab, d) in enumerate(runs):
    R = json.load(open(d + '/reweight_real.json')); p = np.array(R['p']); lo = np.array(R['lo']); hi = np.array(R['hi']); K = len(p); s = np.arange(K) / (K - 1.); col = 'C%d' % j
    ax[0].plot(s, p, color=col, lw=2, label=lab); ax[0].fill_between(s, lo, hi, color=col, alpha=0.25)
    for h in R['halves']: ax[0].plot(s, h, color=col, lw=0.7, ls=':')
    G = lambda q: -np.log(np.maximum(q, 1e-6)) + np.log(p.max())
    ax[1].plot(s, G(p), color=col, lw=2, label=lab); ax[1].fill_between(s, G(hi), G(lo), color=col, alpha=0.25)
ax[0].set_xlabel('path position s (0 = closed conformer 18, 1 = open end)'); ax[0].set_ylabel('occupancy per path step'); ax[0].legend(fontsize=8)
ax[1].set_xlabel('path position s'); ax[1].set_ylabel('apparent G (kT, raw: true heights ~1.5-2x larger on truth)'); ax[1].set_ylim(-0.3, 7)
fig.suptitle('RNase P closed family r1 (92.5k): occupancy along closed -> open morph paths (pose-marginalised, lam 0.01)', fontsize=10); fig.tight_layout(); fig.savefig(out, dpi=120); print('wrote', out)
