#!/usr/bin/env python3
"""Stage the supported zoomed sub-state maps of a run into one flat folder with a ranked README.
usage: stage_zoom_maps.py RUNDIR OUTDIR"""
import re, os, sys, glob, shutil, numpy as np
R, S = sys.argv[1:3]; os.makedirs(S, exist_ok=True)
kept = [l.split()[0] for l in open(R + '/merge_tree_zoom.txt') if re.match(r'f\d_', l)]
arm = {}
for l in open(R + '/arm_zoom.txt'):
    if l.startswith('#'): continue
    w = l.split(); m = re.search(r'fam(\d)/zoomregion([^/]+)/level1/states_dir/state_(\d+)', w[0]); arm['f%s_%s_s%02d' % (m.group(1), m.group(2), int(m.group(3)))] = float(w[1])
ids = [int(x) for x in open(R + '/readout55_zoom/matches.txt').readline().split('|')[1].split()[3:]]
names = [l.strip() for l in open(R + '/readout55_zoom/map_order.txt')]; best = {}
for l in open(R + '/readout55_zoom/matches.txt'):
    if l.startswith('#') or not l.strip(): continue
    w = l.split(); C = np.array([float(x) for x in w[3:]]); o = np.argsort(-C); nm = names[int(w[0]) - 1]
    m = re.search(r'state_fam(\d)_zoomregion(.+)_(\d+)\.mrc', nm); key = 'f%s_%s_s%02d' % (m.group(1), m.group(2), int(m.group(3)))
    best[key] = (ids[o[0]], C[o[0]], C[o[0]] - C[o[1]])
rows = []
for k in kept:
    m = re.match(r'f(\d)_(.+)_s(\d+)$', k); fam, reg, s = m.group(1), m.group(2), int(m.group(3)); d = '%s/fam%s/zoomregion%s' % (R, fam, reg); n = 0
    for l in open(d + '/level1/clusters.txt'):
        if not l.startswith('#') and int(l.split()[0]) == s: n = int(l.split()[3])
    out = 'zoom_f%s_%s_sub%02d.mrc' % (fam, reg, s); shutil.copy('%s/level1/states_dir/state_%02d.mrc' % (d, s), S + '/' + out); b = best.get(k, (0, 0.0, 0.0)); rows.append((out, n, b[0], b[1], b[2], arm.get(k, float('nan'))))
rows.sort(key=lambda r: (-int(bool(r[3] >= 0.4 and r[4] >= 0.05)), -float(r[3])))
with open(S + '/README.txt', 'w') as o:
    o.write('f692c zoomed sub-states: %d maps that passed the half-map support filter (FSC <= 12 A, >= 1000 particles); raw unsharpened reconstructions, box 160, 1.464 A/px.\n' % len(rows))
    o.write('Name: family, region of the zoomed leaf (_1 = the whole-family region, motion_k = motion region k), the leaf (its level-1 state id), the sub-state.\n')
    o.write('Sorted: conformer-matching states first (deviation corr >= 0.40 with a 0.05 margin over every other of the 55 docked apo maps), then by corr. closed_occ: 1.00 = deposited closed map, open conformers 0.11-0.34.\n\n')
    o.write('%-36s %9s %9s %6s %6s %10s\n' % ('map', 'particles', 'best_conf', 'corr', 'margin', 'closed_occ'))
    for r in rows: o.write('%-36s %9d %9d %6.2f %6.2f %10.2f\n' % r)
print('%d maps; %d conformer-matching; %d with closed occupancy >= 0.85' % (len(rows), sum(1 for r in rows if bool(r[3] >= 0.4 and r[4] >= 0.05)), sum(1 for r in rows if r[5] >= 0.85)))
