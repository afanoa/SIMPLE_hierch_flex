#!/usr/bin/env python3
"""Where do the particles of earlier pure closed sub-states go in a new real-data tree? proxy = meta rows carrying the
given labels in an earlier labelling (f692c family-1 zoom sub-states 1 and 3, closed-arm occupancy 0.92); per node and
half: the proxy share of the node and of each child; with FINAL (a final-pass dir) also per final state.
usage: closed_proxy_track.py LABELS_META LABELS(1:3) RUN_PREFIX TAG SIM_PREFIX [FINAL_DIR]"""
import sys, os, glob, numpy as np
labf, labs, runp, tag, simp = sys.argv[1:6]; fin = sys.argv[6] if len(sys.argv) > 6 else ''; want = [int(v) for v in labs.split(':')]
L = np.loadtxt(labf).astype(int); proxy = np.isin(L, want)
for H in 'AB':
    mr = np.loadtxt('%s_%s/meta_rows.txt' % (simp, H)).astype(int) - 1; px = proxy[mr]; base = '%s_%s/%s' % (runp, H, tag)
    nodes = sorted([os.path.basename(d) for d in glob.glob(base + '/r*') if os.path.isfile(d + '/subset.txt')], key=lambda s: [int(x) for x in s.split('.')[1:]])
    print('half %s: %d proxy-closed particles of %d (%.2f %%)' % (H, px.sum(), len(px), 100 * px.mean()), flush=True)
    for nd in nodes:
        sub = np.loadtxt(base + '/' + nd + '/subset.txt', dtype=int, ndmin=1) - 1; s = '  %-10s n %6d proxy %5d (%.2f %%)' % (nd, len(sub), px[sub].sum(), 100 * px[sub].mean())
        ch = base + '/' + nd + '/children.txt'
        if os.path.isfile(ch):
            c = np.loadtxt(ch).astype(int)
            if c.max() > 1: s += ' -> children ' + ' '.join('%d:%d/%.1f%%' % (k, px[c == k].sum(), 100 * px[c == k].mean()) for k in range(1, c.max() + 1))
        if os.path.isfile(base + '/' + nd + '/junk.flag'): s += ' [junk leaf]'
        print(s)
    if fin and os.path.isfile('%s/children_%s.txt' % (fin, H)):
        c = np.loadtxt('%s/children_%s.txt' % (fin, H)).astype(int)
        print('  final states: ' + ' '.join('%d:%d of %d (%.1f%%)' % (k, px[c == k].sum(), (c == k).sum(), 100 * px[c == k].mean()) for k in range(1, c.max() + 1)))
