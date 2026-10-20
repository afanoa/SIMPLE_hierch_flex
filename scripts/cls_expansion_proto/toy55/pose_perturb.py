#!/usr/bin/env python3
"""Pose error of a refined table against the simulation truth (non-junk: direction, in-plane, shift), and truth tables
perturbed by a given direction/in-plane sigma (deg) for the likelihood-ceiling tolerance curve.
usage: pose_perturb.py TRUE REFINED GT_PARTICLES OUTPREFIX SIG1 [SIG2 ...]"""
import sys, numpy as np
truef, reff, gtf, outp = sys.argv[1:5]; sigs = [float(v) for v in sys.argv[5:]]
def oris(fn, keys):
    return np.array([[float(dict(t.split('=', 1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn)])
K = ('e1', 'e2', 'e3', 'x', 'y'); T = oris(truef, K); Rr = oris(reff, K); conf = np.loadtxt(gtf, comments='#')[:, 1].astype(int); ok = conf > 0
def nvec(e1, e2): a, b = np.radians(e1), np.radians(e2); return np.c_[np.cos(a) * np.sin(b), np.sin(a) * np.sin(b), np.cos(b)]
ang = np.degrees(np.arccos(np.clip(np.sum(nvec(T[:, 0], T[:, 1]) * nvec(Rr[:, 0], Rr[:, 1]), 1), -1, 1)))
de3 = np.abs((T[:, 2] - Rr[:, 2] + 180) % 360 - 180); dsh = np.hypot(T[:, 3] - Rr[:, 3], T[:, 4] - Rr[:, 4])
print('refined vs truth, non-junk (%d): direction error median %.2f mean %.2f deg (90th pct %.1f); in-plane median %.2f deg; shift median %.2f px (90th %.1f)' % (ok.sum(), np.median(ang[ok]), ang[ok].mean(), np.percentile(ang[ok], 90), np.median(de3[ok]), np.median(dsh[ok]), np.percentile(dsh[ok], 90)))
rng = np.random.default_rng(1); lines = open(truef).read().split('\n')
for s in sigs:
    out = []
    for i, l in enumerate(lines):
        if not l.strip(): continue
        kv = dict(t.split('=', 1) for t in l.split() if '=' in t)
        e2 = float(kv['e2']) + rng.normal(0, s); e1 = float(kv['e1']) + rng.normal(0, s) / max(abs(np.sin(np.radians(float(kv['e2'])))), 0.2); e3 = float(kv['e3']) + rng.normal(0, s)
        kv['e1'], kv['e2'], kv['e3'] = '%.4f' % (e1 % 360), '%.4f' % np.clip(e2, 0, 180), '%.4f' % (e3 % 360)
        out.append(' '.join('%s=%s' % (k, v) for k, v in kv.items()))
    open('%s_sig%g.txt' % (outp, s), 'w').write('\n'.join(out) + '\n'); print('wrote %s_sig%g.txt' % (outp, s))
