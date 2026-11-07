#!/usr/bin/env python3
"""Shape of the refined-pose error on a synthetic set: percentiles and tail fractions of the direction error (non-junk),
and the error of the perturbed START poses for comparison. usage: pose_err_hist.py TRUE REFINED START GT_PARTICLES"""
import sys, numpy as np
truef, reff, startf, gtf = sys.argv[1:5]
def oris(fn, keys): return np.array([[float(dict(t.split('=', 1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn)])
def nvec(e1, e2): a, b = np.radians(e1), np.radians(e2); return np.c_[np.cos(a) * np.sin(b), np.sin(a) * np.sin(b), np.cos(b)]
T = oris(truef, ('e1', 'e2')); conf = np.loadtxt(gtf, comments='#')[:, 1].astype(int); ok = conf > 0
for name, fn in (('refined', reff), ('start', startf)):
    R = oris(fn, ('e1', 'e2')); n = min(len(R), len(T)); ang = np.degrees(np.arccos(np.clip(np.sum(nvec(T[:n, 0], T[:n, 1]) * nvec(R[:n, 0], R[:n, 1]), 1), -1, 1)))[ok[:n]]
    print('%-8s direction error pct 50/75/90/95/99: %s; fraction > 10 deg %.3f, > 20 %.3f, > 60 %.3f, > 120 %.3f' % (name, ' '.join('%.1f' % v for v in np.percentile(ang, (50, 75, 90, 95, 99))), (ang > 10).mean(), (ang > 20).mean(), (ang > 60).mean(), (ang > 120).mean()))
