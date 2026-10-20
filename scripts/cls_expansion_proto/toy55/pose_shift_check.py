#!/usr/bin/env python3
"""How far did f692e's per-family refine3D_states move the poses away from the consensus? For the family members:
rotation angle between the consensus and refined orientation, the viewing-direction change, and the shift change;
overall and for the closed zoom sub-state particles. usage: pose_shift_check.py META F692E_DIR ZOOM_STATES"""
import sys, numpy as np
meta, E, zf = sys.argv[1:4]
M = np.load(meta); Z = np.loadtxt(zf).astype(int)
def rot(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, s1, c2, s2, c3, s3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -s1; Z1[:, 1, 0] = s1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = s2; Y2[:, 2, 0] = -s2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -s3; Z3[:, 1, 0] = s3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
for fam in (1, 2):
    sub = np.loadtxt('%s/expand/family%d_subset.txt' % (E, fam), dtype=int) - 1; ref = {}
    keys = ('e1', 'e2', 'e3', 'x', 'y'); R = np.zeros((len(M['e1']), 5)); ok = np.zeros(len(M['e1']), bool)
    for i, l in enumerate(open('%s/fam%d/refined_ptcl3D.txt' % (E, fam))):
        if i in ref: pass
    # parse only the members (the file has every meta row)
    want = np.zeros(len(M['e1']), bool); want[sub] = True
    for i, l in enumerate(open('%s/fam%d/refined_ptcl3D.txt' % (E, fam))):
        if want[i]: d = kv(l); R[i] = [float(d.get(k, 0)) for k in keys]; ok[i] = True
    s = sub[ok[sub]]; Rc = rot(M['e1'][s], M['e2'][s], M['e3'][s]); Rr = rot(R[s, 0], R[s, 1], R[s, 2])
    tr = np.einsum('nij,nij->n', Rc, Rr); ang = np.degrees(np.arccos(np.clip((tr - 1) / 2, -1, 1)))
    vc = Rc[:, 2, :]; vr = Rr[:, 2, :]; dirang = np.degrees(np.arccos(np.clip(np.sum(vc * vr, 1), -1, 1)))
    dsh = np.hypot(R[s, 3] - M['x'][s], R[s, 4] - M['y'][s])
    def pc(v): return ' '.join('%.1f' % x for x in np.percentile(v, (25, 50, 75, 90)))
    print('family %d: %d members; rotation change deg p25/50/75/90 %s; viewing-direction change %s; shift change px %s' % (fam, len(s), pc(ang), pc(dirang), pc(dsh)))
    for nm, sel in (('closed pure (zoom 1)', Z[s] == 1), ('closed near (zoom 143/77/30)', np.isin(Z[s], [143, 77, 30]))):
        if sel.any(): print('   %s: %d particles; rotation change %s; viewing-direction change %s; shift change %s' % (nm, sel.sum(), pc(ang[sel]), pc(dirang[sel]), pc(dsh[sel])))
