#!/usr/bin/env python3
"""qtoyflex_maps.py: the maps and conformer table of the flex <-> refinement positive control. qtoy2's six arm conformers (ordered by
s1) become a dense continuum: NLEV maps evenly spaced in s1 over their range, each the linear interpolation of the two bracketing
conformers, with the image weight interpolated the same way (total = qtoy2's arm total). The closed (18) and holo (101) states stay
discrete, as in qtoy2, so the truth is three discrete states (arm, closed, holo) plus junk, with s1 continuous inside the arm.
usage: qtoyflex_maps.py QTOY2_CSV MAPS160 OUTSUB OUT_CSV [NLEV]"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, write_mrc

src, M, sub, out = sys.argv[1:5]
nlev = int(sys.argv[5]) if len(sys.argv) > 5 else 30
rows = list(csv.DictReader(open(src)))
arm = sorted([r for r in rows if r['kind'] == 'arm'], key=lambda r: float(r['s1']))
s = np.array([float(r['s1']) for r in arm])
w = np.array([float(r['images']) for r in arm])
ca = np.array([float(r['closed_arm']) for r in arm])
vols = [read_mrc(os.path.join(M, r['map'])) for r in arm]
os.makedirs(os.path.join(M, sub), exist_ok=True)
lev = np.linspace(s[0], s[-1], nlev)
wl = np.interp(lev, s, w)
wl *= w.sum() / wl.sum()
table = []
for j, t in enumerate(lev):
    k = min(int(np.searchsorted(s, t, side='right')) - 1, len(s) - 2)
    f = (t - s[k]) / (s[k + 1] - s[k])
    name = os.path.join(sub, 'arm%03d_m160.mrc' % (201 + j))
    write_mrc(os.path.join(M, name), ((1. - f) * vols[k] + f * vols[k + 1]).astype(np.float32), 1.464)
    table.append(dict(conformer=201 + j, images='%.4f' % wl[j], kind='arm', s1='%.2f' % t,
                      closed_arm='%.3f' % np.interp(t, s, ca), map=name))
table += [{k: r[k] for k in ('conformer', 'images', 'kind', 's1', 'closed_arm', 'map')} for r in rows if r['kind'] != 'arm']
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=['conformer', 'images', 'kind', 's1', 'closed_arm', 'map'], lineterminator='\n')
    wr.writeheader()
    wr.writerows(table)
print('%d arm levels, s1 %.2f..%.2f, weights %.2f..%.2f; discrete: %s' % (nlev, lev[0], lev[-1], wl.min(), wl.max(),
      ', '.join('%s %s (%s)' % (r['kind'], r['conformer'], r['images']) for r in rows if r['kind'] != 'arm')))
