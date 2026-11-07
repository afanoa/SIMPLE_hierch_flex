#!/usr/bin/env python3
"""Synthetic pose accuracy: for the rows a labelling covers (LABELS > 0), the rotation error of every pose table against the
simulation's truth (all.txt), after removing one global frame rotation per table; viewing-direction error and shift error;
overall, open vs closed conformers and conformer 18, and how often a table beats the first one.
usage: pose_truth_check.py TRUTH GT LABELS NAME=TABLE [NAME=TABLE ...]"""
import sys, re, numpy as np
truth, gtf, labf = sys.argv[1:4]; tabs = [t.split('=', 1) for t in sys.argv[4:]]
lab = np.loadtxt(labf).astype(int); rows = np.where(lab > 0)[0]; want = np.zeros(len(lab), bool); want[rows] = True
gt = np.loadtxt(gtf, comments='#').astype(float); conf = gt[:, 1].astype(int); closed = gt[:, 2].astype(int)
pat = {k: re.compile(r'(?:^|\s)%s=(\S+)' % k) for k in ('e1', 'e2', 'e3', 'x', 'y')}
def read(f):
    out = np.zeros((len(lab), 5)); i = 0
    for l in open(f):
        if not l.strip(): continue
        if i < len(lab) and want[i]: out[i] = [float(m.group(1)) if (m := pat[k].search(l)) else 0. for k in ('e1', 'e2', 'e3', 'x', 'y')]
        i += 1
    return out[rows]
def rotm(e1, e2, e3):   # splitbench5.rotm (calibrated convention: p3 = R k)
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, s1, c2, s2, c3, s3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -s1; Z1[:, 1, 0] = s1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = s2; Y2[:, 2, 0] = -s2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -s3; Z3[:, 1, 0] = s3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
Tp = read(truth); Rt = rotm(Tp[:, 0], Tp[:, 1], Tp[:, 2]); c = conf[rows]; cl = closed[rows]
groups = (('all', np.ones(len(rows), bool)), ('open', (c > 0) & (cl == 0)), ('closed', cl == 1), ('conf18', c == 18), ('junk', c == 0))
pc = lambda v: '%.2f/%.2f/%.2f' % tuple(np.percentile(v, (50, 75, 90)))
print('%d rows (labels > 0); errors median/p75/p90' % len(rows)); first = None
for nm, f in tabs:
    P = read(f); Re = rotm(P[:, 0], P[:, 1], P[:, 2]); U, _, Vt = np.linalg.svd(np.einsum('nij,nkj->ik', Re, Rt)); G = U @ np.diag([1, 1, np.linalg.det(U @ Vt)]) @ Vt
    g = np.degrees(np.arccos(np.clip((np.trace(G) - 1) / 2, -1, 1))); RtG = np.einsum('ij,njk->nik', G, Rt)
    ang = np.degrees(np.arccos(np.clip((np.einsum('nij,nij->n', Re, RtG) - 1) / 2, -1, 1)))
    vd = np.degrees(np.arccos(np.clip(np.sum(Re[:, :, 2] * RtG[:, :, 2], 1), -1, 1)))
    sh = min((np.hypot(P[:, 3] - s * Tp[:, 3], P[:, 4] - s * Tp[:, 4]) for s in (1, -1)), key=np.median)
    print('%-14s frame offset %.2f deg' % (nm, g))
    for gn, m in groups:
        if m.sum() < 50: continue
        extra = '' if first is None else '; better than %s %.0f %%, worse %.0f %%' % (first[0], 100 * np.mean(ang[m] < first[1][m] - 0.5), 100 * np.mean(ang[m] > first[1][m] + 0.5))
        print('   %-7s n %6d  rotation %s  view dir %s  shift px %s%s' % (gn, m.sum(), pc(ang[m]), pc(vd[m]), pc(sh[m]), extra))
    if first is None: first = (nm, ang)
