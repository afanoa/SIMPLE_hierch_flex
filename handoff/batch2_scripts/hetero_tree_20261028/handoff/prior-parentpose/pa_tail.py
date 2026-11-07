#!/usr/bin/env python3
# where the large-jump tail of the child's pose drift (incoming -> stage B) falls by truth class.
# usage: pa_tail.py LABEL NODE TRUTH N_DUMP B_DUMP
import sys, os, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent')
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor')
import latent_states_v4 as LS
import importlib.util
spec = importlib.util.spec_from_file_location('pd', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/pa_drift_lib.py')
lab, node, truth, fn, fb = sys.argv[1:6]
exec(open('/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/pa_drift_lib.py').read())
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
kind, T = LS.load_truth(truth, node, rows)
N, B = read(fn), read(fb)
m = np.ones(rows.size, bool)
geo, ddir, tw, dsh = change(N, B, m, 1.464)
tail = (tw > 30) | (dsh > 10)
st = B['state'].astype(int)
print('## %s: n %d, tail (in-plane > 30 deg or shift > 10 A) %d = %.4f; direction > 8 deg %.4f' % (lab, rows.size, tail.sum(), tail.mean(), (ddir > 8).mean()))
if kind == 'closed':
    groups = {'closedref': T['closedref'] > 0.5, 'other': T['closedref'] <= 0.5}
else:
    j = T['junk'] > 0
    groups = {'junk': j, 'closed(18)': (T['combo'] == 18) & ~j, 'holo(101)': (T['combo'] == 101) & ~j}
    groups['rest'] = ~(j | groups['closed(18)'] | groups['holo(101)'])
for g, mk in groups.items():
    if mk.sum() == 0:
        continue
    print('  %-11s n %6d (%.3f of node)  tail %.4f  share of tail %.3f  | in-plane p99 %.1f  shift p99 %.1f A  direction med %.2f p90 %.2f' %
          (g, mk.sum(), mk.mean(), tail[mk].mean(), (tail & mk).sum() / max(1, tail.sum()), np.percentile(tw[mk], 99), np.percentile(dsh[mk], 99),
           np.median(ddir[mk]), np.percentile(ddir[mk], 90)))
print('  stage-B state of the tail: %s; of all: %s' % (' '.join('%d:%d' % (u, (st[tail] == u).sum()) for u in np.unique(st)),
                                                      ' '.join('%d:%d' % (u, (st == u).sum()) for u in np.unique(st))))
