#!/usr/bin/env python3
"""Arm-state readout of maps: closed-arm mask = where the deposited closed conformer (18) has density the consensus lacks,
open-arm mask = the reverse (both from the 10 A low-passed, standardised maps); per map the occupancy of each mask relative
to its reference (closed map / consensus) -> closed, open or intermediate (partial occupancy of both = swinging/mixed).
usage: arm_readout.py CONSENSUS CLOSEDMAP OUT.txt MAP [MAP ...]  (MAP may be a directory: all state_*.mrc / recvol_state??.mrc)"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
cons_f, closed_f, out = sys.argv[1:4]; inputs = sys.argv[4:]; smpd = 1.464; LP = 10.
n = read_mrc(cons_f).shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G))
cons = lp(read_mrc(cons_f)); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(read_mrc(closed_f))); d = closed - cons; thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99)
mC = d > thr; mO = d < -thr; refC = closed[mC].mean(); refO = cons[mO].mean()
print('closed-arm mask %d voxels (%.1f%% of the support), open-arm mask %d voxels; reference occupancies %.2f / %.2f' % (mC.sum(), 100. * mC.sum() / supp.sum(), mO.sum(), refC, refO))
files = []
for p in inputs:
    files += sorted(glob.glob(p + '/state_*.mrc')) + sorted(glob.glob(p + '/recvol_state??.mrc')) if os.path.isdir(p) else [p]
rows = []
for fn in files:
    v = std(lp(read_mrc(fn))); c = v[mC].mean() / refC; o = v[mO].mean() / refO
    lab = 'closed' if (c > 0.6 and o < 0.45) else 'open' if (o > 0.6 and c < 0.45) else 'intermediate'
    rows.append((fn, c, o, lab)); print('%-70s closed %.2f open %.2f -> %s' % (fn.split('rnasep_klem/')[-1], c, o, lab))
with open(out, 'w') as f:
    f.write('# map closed_arm_occupancy open_arm_occupancy call (occupancies relative to the deposited closed map / the consensus)\n')
    for r in rows: f.write('%s %.3f %.3f %s\n' % r)
calls = [r[3] for r in rows]; print('summary: %d closed, %d open, %d intermediate of %d maps' % (calls.count('closed'), calls.count('open'), calls.count('intermediate'), len(calls)))
