#!/usr/bin/env python3
"""land_plane.py: the arm-readout plane (closed-site occupancy vs open-site occupancy, arm_readout.py's definitions). Points: tree leaves (plain
maps; area ~ particles), deposited apo / holo conformers, rigid-swing maps (a hinge rotation of the closed arm), optional window maps along a
coordinate (connected in order). Writes PNG + the particle-weighted marginal of the leaves' closed-site occupancy.
usage: land_plane.py OUT.png --leaves ARM.txt --counts UNITS_LIST --confs ARM_CONFS.txt [--windows ARM.txt:label ...]"""
import argparse, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--leaves', required=True); ap.add_argument('--counts', required=True); ap.add_argument('--confs', required=True)
ap.add_argument('--windows', nargs='*', default=[]); a = ap.parse_args()
def tab(f): return [(os.path.basename(l.split()[0]), float(l.split()[1]), float(l.split()[2])) for l in open(f) if l.strip() and not l.startswith('#')]
L = tab(a.leaves); n = np.array([float(l.split()[3]) for l in open(a.counts) if l.strip() and not l.startswith('#')]); C = tab(a.confs)
fig, ax = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
x = ax[0]; lc = np.array([l[1] for l in L]); lo = np.array([l[2] for l in L])
x.scatter(lo, lc, s=n / n.max() * 400, c='C0', alpha=0.45, edgecolor='k', lw=0.4, label='tree leaves (plain maps, area ~ particles)')
def cid(nm): return int(nm[4:].split('_')[0])
for kind, lo_, hi_, col, mk in (('deposited apo', 0, 99, 'C3', 'x'), ('deposited holo', 100, 199, 'C1', '+')):
    pts = np.array([(c[2], c[1]) for c in C if '_emd' in c[0] and lo_ <= cid(c[0]) <= hi_]); x.scatter(pts[:, 0], pts[:, 1], marker=mk, c=col, s=28, label=kind)
for lab, lo_, hi_, col in (('rigid swing 0-30 A', 200, 230, 'k'), ('rigid swing 32-60 A', 400, 430, '0.6')):
    sw = np.array(sorted([(cid(c[0]), c[2], c[1]) for c in C if '_swing' in c[0] and lo_ <= cid(c[0]) <= hi_]))
    if len(sw): x.plot(sw[:, 1], sw[:, 2], '-', color=col, lw=1.2, label=lab)
for w in a.windows:
    f, lab = w.split(':'); W = tab(f); x.plot([q[2] for q in W], [q[1] for q in W], 'o-', ms=4, lw=1.4, label=lab)
x.set_xlabel('open-site occupancy (rel. consensus)'); x.set_ylabel('closed-site occupancy (rel. deposited closed conf18)'); x.legend(fontsize=7, loc='upper right'); x.grid(alpha=0.3)
x = ax[1]; bins = np.linspace(-0.1, 1.2, 27); h = np.histogram(lc, bins=bins, weights=n)[0] / n.sum()
x.bar(0.5 * (bins[1:] + bins[:-1]), h, width=np.diff(bins), color='C0', alpha=0.6, edgecolor='k', lw=0.4)
x.set_xlabel('closed-site occupancy of the leaf map'); x.set_ylabel('fraction of particles (leaf granularity)')
x.set_title('map-level marginal (leaves are mixtures: not deconvolved)', fontsize=9)
fig.tight_layout(); fig.savefig(a.out, dpi=120)
print('leaves %d, particles %d; particle-weighted closed-site occupancy: mean %.3f; share of particles in leaves >= 0.8: %.3f, 0.45-0.8: %.3f, < 0.45: %.3f' % (
    len(L), n.sum(), (lc * n).sum() / n.sum(), n[lc >= 0.8].sum() / n.sum(), n[(lc >= 0.45) & (lc < 0.8)].sum() / n.sum(), n[lc < 0.45].sum() / n.sum()))
