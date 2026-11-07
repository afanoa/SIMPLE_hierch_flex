#!/usr/bin/env python3
# which particles an arm relabels relative to production: by truth class, by production's last-iteration state stability
# (mi_state = 0: the particle changed state in production's last iteration) and by the production drift tail.
# usage: pa_change.py KEY NODE TRUTH ARM_B_DUMP [ARM_B_DUMP ...]
import sys, os, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent')
import latent_states_v4 as LS
exec(open('/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/pa_drift_lib.py').read())
A = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor'
key, node, truth = sys.argv[1:4]
rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
kind, T = LS.load_truth(truth, node, rows)
def mi(fn):
    v = []
    for line in open(fn):
        d = dict(t.split('=', 1) for t in line.split() if '=' in t)
        if 'e1' in d: v.append(float(d.get('mi_state', 'nan')))
    return np.array(v)
P = read(os.path.join(A, 'drift', key + '_b.txt')); N = read(os.path.join(A, 'drift', key + '_n.txt'))
mis = mi(os.path.join(A, 'drift', key + '_b.txt'))
geo, ddir, tw, dsh = change(N, P, np.ones(rows.size, bool), 1.464)
tail = (tw > 30) | (dsh > 10)
cls = (T['closedref'] > 0.5) if kind == 'closed' else (T['junk'] > 0)
cname = 'closedref' if kind == 'closed' else 'junk'
print('## %s: production last-iteration state flippers (mi_state=0) %.3f, tail %.3f, %s %.3f' % (key, (mis == 0).mean(), tail.mean(), cname, cls.mean()))
for fb in sys.argv[4:]:
    s = read(fb)['state'].astype(int); p = P['state'].astype(int)
    if len(np.unique(s)) == 2 and len(np.unique(p)) == 2 and (s != p).mean() > 0.5: s = 3 - s
    ch = s != p
    print('  %s: relabelled %.3f | among flippers %.3f, non-flippers %.3f | tail %.3f, non-tail %.3f | %s %.3f, rest %.3f | share of relabelled that were flippers %.2f' %
          (fb.split('/')[-3], ch.mean(), ch[mis == 0].mean(), ch[mis == 1].mean(), ch[tail].mean() if tail.any() else np.nan, ch[~tail].mean(),
           cname, ch[cls].mean(), ch[~cls].mean(), (ch & (mis == 0)).sum() / max(1, ch.sum())))
