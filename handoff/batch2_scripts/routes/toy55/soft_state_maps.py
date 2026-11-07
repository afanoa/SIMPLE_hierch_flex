#!/usr/bin/env python3
"""Soft state maps and per-particle weights of a board's states (no reconstruction): a state = a cluster of cells
(clusters.txt); its soft map = consensus + the mean deviation of its member cells (klem_halves.npz D0/D1 are the
cross-half cell deviations on the board's Fourier cube), written per half and merged on the full box; its weight for
every particle = the summed responsibility over its member cells (klem_resp.npz). The Fourier placement convention is
checked against a written cell map (cells/cell_*.mrc) and the better of the two shift conventions is used.
usage: soft_state_maps.py BOARD CLUSTERS.txt CONSENSUS OUTDIR [--smpd 1.464]"""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('clusters'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); H = np.load(a.board + '/klem_halves.npz'); R = np.load(a.board + '/klem_resp.npz'); D0, D1, band, n3 = H['D0'], H['D1'], H['band'], int(H['n3']); r, pinds = R['r'], R['pinds'].astype(int)
cons = read_mrc(a.cons).astype(np.float32); N = cons.shape[0]; c = N // 2; h = n3 // 2
kk_ = np.arange(n3) - n3 // 2; chk = ((-1.0) ** (kk_[:, None, None] + kk_[None, :, None] + kk_[None, None, :])).astype(np.float32)
Fc = np.fft.fftshift(np.fft.fftn(cons)); Fc_crop = Fc[c - h:c + h + 1, c - h:c + h + 1, c - h:c + h + 1]
def place(dev, conv):
    cube = np.zeros((n3, n3, n3), np.complex64); cube.ravel()[band] = dev
    if conv & 1: cube = cube[::-1, ::-1, ::-1]
    if conv & 2: cube = cube.transpose(2, 1, 0)
    if conv & 4: cube = np.conj(cube)
    if conv & 8: cube = cube * chk   # phase origin at the box centre: (-1)^(h+k+l)
    return cube
def tomap(dev, conv, scale):
    F = Fc.copy(); F[c - h:c + h + 1, c - h:c + h + 1, c - h:c + h + 1] = Fc_crop + place(dev, conv) * scale
    return np.real(np.fft.ifftn(np.fft.ifftshift(F))).astype(np.float32)
# convention check on the DEVIATION of a written cell map (cell map - consensus) against the placed D of that cell
cellfiles = sorted(glob.glob(a.board + '/cells/cell_*.mrc')); conv, bscale = 0, 1.0
cellmaps = [read_mrc(f).astype(np.float32) for f in cellfiles]
if len(cellfiles) >= 2:   # convention from a cell DIFFERENCE (the board's own consensus cancels): F(cell_1 - cell_2) vs placed (D_1 - D_2), in the band
    Fr = np.fft.fftshift(np.fft.fftn(cellmaps[0] - cellmaps[1]))[c - h:c + h + 1, c - h:c + h + 1, c - h:c + h + 1]; dd = 0.5 * (D0[0] + D1[0]) - 0.5 * (D0[1] + D1[1]); best = -2
    for cv in range(16):
        P = place(dd, cv); cc = np.vdot(P, Fr).real / max(np.sqrt(np.vdot(P, P).real * np.vdot(Fr, Fr).real), 1e-30)
        if cc > best: best, conv = cc, cv
    P = place(dd, conv); bscale = float(np.vdot(P, Fr).real / max(np.vdot(P, P).real, 1e-30))
    print('deviation convention check on cell_1 - cell_2 (Fourier, in-band): corr %.3f (conv %d, scale %.4g)%s' % (best, conv, bscale, '' if best > 0.95 else '  WARNING: no convention reproduces the cell difference'))
def half_from_full(full, dhalf, dmean):   # half map = full soft map + placed (D_half - D_mean)
    F = np.fft.fftshift(np.fft.fftn(full)); F[c - h:c + h + 1, c - h:c + h + 1, c - h:c + h + 1] += place(dhalf - dmean, conv) * bscale
    return np.real(np.fft.ifftn(np.fft.ifftshift(F))).astype(np.float32)
states = {}
for l in open(a.clusters):
    if l.startswith('#'): continue
    p = l.split(); states[int(p[0])] = np.array([int(x) - 1 for x in p[4:]])
for s, cells in sorted(states.items()):
    full = np.mean([cellmaps[k] for k in cells], 0).astype(np.float32); d0, d1 = D0[cells].mean(0), D1[cells].mean(0); dm = 0.5 * (d0 + d1)
    write_mrc('%s/state_soft_%02d.mrc' % (a.out, s), full, a.smpd); write_mrc('%s/state_soft_%02d_even.mrc' % (a.out, s), half_from_full(full, d0, dm), a.smpd); write_mrc('%s/state_soft_%02d_odd.mrc' % (a.out, s), half_from_full(full, d1, dm), a.smpd)
    w = r[:, cells].sum(1); np.savetxt('%s/weights_%02d.txt' % (a.out, s), np.c_[pinds, w], fmt='%d %.5f')
    print('state %d: %d cells, weight mass %.0f, %d particles with weight > 0.5 -> state_soft_%02d.mrc (+even/odd), weights_%02d.txt' % (s, len(cells), w.sum(), (w > 0.5).sum(), s, s))
