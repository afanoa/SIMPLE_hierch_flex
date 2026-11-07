#!/usr/bin/env python3
"""lsamp_poleaxis.py: a run's landscape coordinate taken as the POLE AXIS instead of the leading score axis: per particle the score
against the closed-pole map minus the score against the open-pole map (the log-likelihood ratio of the two poles at the run's poses,
per pixel set, divided by the particle's noise level as lsamp2h.py --ppnorm 1), standardised and written as coords.npz (tA from pixel set A,
tB from set B, tfull from both; closed end low) for lsamp_repro_bins2.py, with an orientation stub arm_final_all.txt (no maps behind it).
usage: lsamp_poleaxis.py --scores SC --closed NAME --open NAME --out OUT"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--scores', required=True); ap.add_argument('--closed', required=True); ap.add_argument('--open', required=True); ap.add_argument('--out', required=True)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
M = np.load(a.scores + '/meta.npz'); names = [str(x) for x in M['names']]; EO = M['eo']; OWN = M['own']; rows = np.where((OWN >= 0) & (EO >= 0))[0]; ic, io = names.index(a.closed), names.index(a.open)
SA = np.load(a.scores + '/scores_pA.npy', mmap_mode='r')[rows]; SB = np.load(a.scores + '/scores_pB.npy', mmap_mode='r')[rows]
ok = np.isfinite(SA).all(1) & np.isfinite(SB).all(1); rows, SA, SB = rows[ok], np.asarray(SA[ok], np.float64), np.asarray(SB[ok], np.float64)
lam = np.asarray(M['ppnoise'])[rows].astype(np.float64); lam = np.where(np.isfinite(lam) & (lam > 0), lam, np.nanmedian(lam))
dA = (SA[:, ic] - SA[:, io]) / lam; dB = (SB[:, ic] - SB[:, io]) / lam
def z(x): return ((x - x.mean()) / x.std())[:, None]
np.savez(os.path.join(a.out, 'coords.npz'), rows=rows, tfull=-z(dA + dB), tA=-z(dA), tB=-z(dB), ax=np.array([0]), lam=np.array([np.nan]), nper=np.array([0]))
open(os.path.join(a.out, 'arm_final_all.txt'), 'w').write('# orientation stub for lsamp_repro_bins2.py (no maps): the coordinate is stored closed end low\nstate_01.mrc 1 0 stub\nstate_02.mrc 0 1 stub\n')
print('%s: %d particles; pole axis %s minus %s; fold A vs fold B coordinate corr %.3f' % (a.scores, len(rows), a.closed, a.open, np.corrcoef(dA, dB)[0, 1]))
