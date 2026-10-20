#!/usr/bin/env python3
"""Truth-initialised lattice positions for a focused region board on the hinged-arm synthetic set: the arm is the one
whose travel (|map(conformer with max theta_k) - map(min theta_k)|) overlaps the region mask most; every particle's
cell is (theta quantile bin of n1, phi quantile bin of n2) of that arm's TRUE angles (junk: random cell). Writes
OUT/init_positions.txt ("pind g1 g2") and OUT/init_labels.txt (one label per particle, 0 = junk) for the readout.
usage: oracle_arm_init.py GT_PARTICLES TRUTH_ARMS MASK SYNTHMAPS OUT [--n1 7 --n2 3 --arm K]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('gt'); ap.add_argument('truth'); ap.add_argument('mask'); ap.add_argument('maps'); ap.add_argument('out')
ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=3); ap.add_argument('--arm', type=int, default=0); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(a.truth)) if w and not w[0].startswith('#')}
G = np.loadtxt(a.gt, comments='#'); pind = G[:, 0].astype(int); conf = G[:, 1].astype(int); narm = len(next(iter(ang.values()))) // 2
msk = read_mrc(a.mask).astype(np.float32); msk = msk / max(msk.max(), 1e-12)
mp = lambda c: read_mrc('%s/conf%02d_synth.mrc' % (a.maps, c)).astype(np.float32)
arm = a.arm
if arm == 0:
    best = None
    for k in range(1, narm + 1):
        th = {c: v[2 * (k - 1)] for c, v in ang.items()}; cmax = max(th, key=th.get); cmin = min(th, key=th.get)
        travel = np.abs(mp(cmax) - mp(cmin)); ov = float((travel * msk).sum() / max(travel.sum(), 1e-12))
        print('arm %d: travel inside the mask %.3f (conformers %d vs %d, theta %.1f vs %.1f)' % (k, ov, cmax, cmin, th[cmax], th[cmin]))
        if best is None or ov > best[0]: best = (ov, k)
    arm = best[1]; print('arm %d chosen (overlap %.3f)' % (arm, best[0]))
th = np.array([ang[c][2 * (arm - 1)] if c > 0 else np.nan for c in conf]); ph = np.array([ang[c][2 * (arm - 1) + 1] if c > 0 else np.nan for c in conf])
ok = conf > 0; rng = np.random.default_rng(1)
def qbin(v, n):   # quantile bins over the non-junk particles -> 0..n-1
    b = np.zeros(len(v), int); q = np.quantile(v[ok], np.linspace(0, 1, n + 1)[1:-1]); b[ok] = np.searchsorted(q, v[ok]); b[~ok] = rng.integers(0, n, (~ok).sum()); return b
g1, g2 = qbin(th, a.n1), qbin(ph, a.n2)
np.savetxt(a.out + '/init_positions.txt', np.c_[pind, g1, g2], fmt='%d')
lab = np.where(ok, g1 + a.n1 * g2 + 1, 0); np.savetxt(a.out + '/init_labels.txt', lab, fmt='%d')
print('init: %d particles, arm %d, %d x %d cells (theta bins x phi bins), %d junk at random cells -> %s' % (len(pind), arm, a.n1, a.n2, (~ok).sum(), a.out))
