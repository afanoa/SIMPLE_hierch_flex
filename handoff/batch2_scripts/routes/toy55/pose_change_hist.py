#!/usr/bin/env python3
"""Where do large pose changes between two tables happen? rotation-angle and viewing-direction change histograms
overall and per state of a labelling (rows with label > 0). usage: pose_change_hist.py BEFORE AFTER LABELS"""
import sys, re, numpy as np
fb, fa, fl = sys.argv[1:4]; lab = np.loadtxt(fl).astype(int); want = lab > 0
pat = {k: re.compile(r'(?:^|\s)%s=(\S+)' % k) for k in ('e1', 'e2', 'e3')}
def read(f):
    out = np.zeros((len(lab), 3)); i = 0
    for l in open(f):
        if not l.strip(): continue
        if i < len(lab) and want[i]: out[i] = [float(m.group(1)) if (m := pat[k].search(l)) else 0. for k in ('e1', 'e2', 'e3')]
        i += 1
    return out[want]
def rotm(e):
    a1, a2, a3 = np.radians(e[:, 0]), np.radians(e[:, 1]), np.radians(e[:, 2]); c1, s1, c2, s2, c3, s3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -s1; Z1[:, 1, 0] = s1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = s2; Y2[:, 2, 0] = -s2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -s3; Z3[:, 1, 0] = s3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
Rb, Ra = rotm(read(fb)), rotm(read(fa)); L = lab[want]
ang = np.degrees(np.arccos(np.clip((np.einsum('nij,nij->n', Rb, Ra) - 1) / 2, -1, 1)))
vd = np.degrees(np.arccos(np.clip(np.sum(Rb[:, :, 2] * Ra[:, :, 2], 1), -1, 1)))
edges = [0, 2, 5, 10, 20, 45, 90, 135, 181]
def h(v): c = np.histogram(v, edges)[0]; return ' '.join('%d-%d:%.1f%%' % (edges[i], edges[i + 1], 100. * c[i] / len(v)) for i in range(len(c)))
print('%d particles; rotation change: %s' % (len(ang), h(ang))); print('   viewing-direction change: %s' % h(vd))
print('   share with rotation change > 20 deg whose viewing direction moved < 5 deg (in-plane only): %.1f %%' % (100. * np.mean(vd[ang > 20] < 5)))
for k in np.unique(L):
    m = L == k; print('   state %d (n %d): rotation > 20 deg %.1f %%, > 90 deg %.1f %%; median %.1f' % (k, m.sum(), 100 * np.mean(ang[m] > 20), 100 * np.mean(ang[m] > 90), np.median(ang[m])))
