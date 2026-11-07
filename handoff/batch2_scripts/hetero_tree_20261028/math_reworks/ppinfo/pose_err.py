#!/usr/bin/env python3
"""pose_err.py: the pose error of a node's input and stage-B poses against the simulation truth (angle between rotations, degrees;
shift error in pixels), for all particles and for the rare ones. usage: pose_err.py DUMPDIR TRUE_ORITAB TRUTH RARE"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ppinfo_site as S
d, tfn, truth, rare = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
nd = S.read_dump(os.path.join(d, 'n_ptcl3D.txt')); bd = S.read_dump(os.path.join(d, 'b_ptcl3D.txt'))
tr = S.read_oritab(tfn); ind = nd['indstk'].astype(int)
t = np.loadtxt(truth, ndmin=2); conf = np.zeros(int(t[:, 0].max()) + 1, int); conf[t[:, 0].astype(int)] = t[:, 7].astype(int)
israre = conf[np.minimum(ind, len(conf) - 1)] == rare
def ang(ea, eb):
    Ra, Rb = S.euler2m(ea), S.euler2m(eb)
    E = np.einsum('nji,njk->nik', Ra, Rb)
    return np.degrees(np.arccos(np.clip((np.trace(E, axis1=1, axis2=2) - 1) / 2, -1, 1)))
for name, e in (('input', np.c_[nd['e1'], nd['e2'], nd['e3'], nd['x'], nd['y']]), ('stageB', np.c_[bd['e1'], bd['e2'], bd['e3'], bd['x'], bd['y']])):
    a = ang(e[:, :3], tr[ind - 1, :3]); sh = np.hypot(e[:, 3] - tr[ind - 1, 3], e[:, 4] - tr[ind - 1, 4])
    for lab, m in (('all', np.ones(len(a), bool)), ('rare', israre), ('majority', ~israre)):
        print('%-7s %-8s n %6d angle median %5.2f mean %5.2f  >5deg %.3f  >10deg %.3f  >30deg %.3f | shift median %.2f px (box 160)' % (
            name, lab, m.sum(), np.median(a[m]), a[m].mean(), (a[m] > 5).mean(), (a[m] > 10).mean(), (a[m] > 30).mean(), np.median(sh[m])))
