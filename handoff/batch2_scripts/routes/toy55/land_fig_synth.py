#!/usr/bin/env python3
"""land_fig_synth.py: synthetic validation figure: per imposed landscape, the truth, the map-level estimate (window unmixing), the naive
histogram of the calibrated coordinate and the oracle-likelihood argmax / point-coordinate histograms (land_e0 at node poses), as
occupancy and as G = -ln p. usage: land_fig_synth.py OUT.png E1DIR E0CURVES.npz"""
import sys, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out, e1d, e0c = sys.argv[1:4]; C0 = np.load(e0c); shapes = ['flat', 'broad', 'dw2', 'dw4', 'shoulder']; tg = np.arange(31)
def sm3(p): return np.convolve(p, np.ones(3), 'same') / np.convolve(np.ones_like(p), np.ones(3), 'same')
fig, ax = plt.subplots(2, len(shapes), figsize=(3.4 * len(shapes), 6.6), squeeze=False)
for j, s in enumerate(shapes):
    R = json.load(open('%s/%s/e1.json' % (e1d, s))); c = R['curves']; pt = np.array(c['truth']); ref = np.log(sm3(pt).max())
    series = [('truth', pt, 'k', 2.4), ('map-level unmixing', np.array(c['map']), 'C3', 1.6), ('naive coordinate histogram', np.array(c['naive']), 'C2', 1.0)]
    k1 = '%s_tau1__hard' % s
    if k1 in C0.files: series.append(('oracle maps, argmax (node poses)', C0[k1], 'C1', 1.0))
    for i in (0, 1):
        x = ax[i, j]
        for lab, p, col, lw in series:
            x.plot(tg, p if i == 0 else -np.log(np.maximum(sm3(p), 1e-6)) + ref, color=col, lw=lw, label=lab)
        x.set_title(s); x.set_xlabel('swing t (A tip travel)'); x.set_ylabel('occupancy' if i == 0 else 'G (kT)')
        if i == 1: x.set_ylim(-0.5, 6)
ax[0, 0].legend(fontsize=7); fig.tight_layout(); fig.savefig(out, dpi=115); print('wrote', out)
