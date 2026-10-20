#!/usr/bin/env python3
"""Per cell of an arm-focused board: particles (argmax), closed-arm occupancy (closed_arm_cells.txt, cell file order) and
the share of the earlier pure closed sub-state particles (proxy); writes OUT_subset.txt = the particles of the cells whose
occupancy is >= --occ (most closed first, until --nmin particles are reached even if fewer cells pass).
usage: cell_report.py BOARD_DIR PROXY_LABELS 1:3 META_ROWS OUT_PREFIX [--occ 0.6 --nmin 3000]"""
import sys, glob, os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('proxy'); ap.add_argument('want'); ap.add_argument('mr'); ap.add_argument('out')
ap.add_argument('--occ', type=float, default=0.6); ap.add_argument('--nmin', type=int, default=3000); a = ap.parse_args()
cfg = dict(l.split()[:2] for l in open(a.board + '/board.txt') if l.strip()); ng = int(cfg['n1']) * int(cfg['n2'])
fs = sorted(glob.glob(a.board + '/cells/*.mrc')); C = np.loadtxt(a.board + '/closed_arm_cells.txt', ndmin=2); occ = dict(zip(C[:, 0].astype(int), C[:, 1]))
L = np.loadtxt(a.board + '/klem_latent.txt', comments='#', ndmin=2); pind = L[:, 0].astype(int); am = L[:, 3].astype(int)
P = np.isin(np.loadtxt(a.proxy).astype(int), [int(v) for v in a.want.split(':')]); mr = np.loadtxt(a.mr).astype(int) - 1; px = P[mr][pind - 1]
print('%d cell files (%s .. %s), %d grid cells; %d particles, proxy %d (%.2f %%)' % (len(fs), os.path.basename(fs[0]), os.path.basename(fs[-1]), ng, len(pind), px.sum(), 100 * px.mean()))
rows = []
for k in range(1, len(fs) + 1):
    m = am == k; rows.append((k, os.path.basename(fs[k - 1]), int(m.sum()), occ.get(k, np.nan), px[m].sum(), 100. * px[m].mean() if m.any() else 0.))
rows.sort(key=lambda r: -r[3] if r[3] == r[3] else 0)
for r in rows[:10]: print('  cell %2d %-16s n %6d closed-arm %.2f proxy %4d (%.1f %%)' % r)
sel = []; tot = 0
for r in rows:
    if r[2] == 0: continue
    if r[3] >= a.occ or tot < a.nmin: sel.append(r[0]); tot += r[2]
    else: break
m = np.isin(am, sel); np.savetxt(a.out + '_subset.txt', np.sort(pind[m]), fmt='%d')
print('selected cells %s: %d particles, proxy %d (%.1f %%), occupancy %.2f-%.2f' % (sel, m.sum(), px[m].sum(), 100 * px[m].mean(), min(occ.get(k, 0) for k in sel), max(occ.get(k, 0) for k in sel)))
