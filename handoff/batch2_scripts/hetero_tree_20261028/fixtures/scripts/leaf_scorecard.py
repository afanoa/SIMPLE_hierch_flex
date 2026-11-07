#!/usr/bin/env python3
"""closed-state leaf scorecard of finished hetero_tree runs: per leaf (hetero_tree_leaves.txt) the reference closed rows (8,720, the
old full692 tree's closed leaf) it holds (count, recall, purity, F1) and the closed-arm / open-arm readout of its plain leaf map
(make_view 10 A views at box 160 + arm_readout; closed corner >= 0.95 / <= 0.6). Per tree: the best closed-leaf F1, the leaves in the
closed corner, and the best union of purity-ranked leaves (does the tree split closed over several leaves?).
usage: leaf_scorecard.py OUT_DIR NAME=TREE_DIR [NAME=TREE_DIR ..]"""
import sys, os, subprocess, numpy as np
from pathlib import Path
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
TREES = [a.split('=', 1) for a in sys.argv[2:]]
SC = Path(__file__).resolve().parent; PY = sys.executable; NR = 932496
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160')
CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
ref = np.zeros(NR + 1, bool)
ref[np.loadtxt('/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt', dtype=int)] = True
spec, leaves = [], []
for name, T in TREES:
    fx, rmap = Path(T).parent / 'fixture_rows.txt', None
    if fx.exists():   # a row-local fixture run: its row j is project row fixture_rows.txt(j)
        rmap = np.loadtxt(fx, dtype=int, ndmin=1)
        if sum(1 for l in open(Path(T) / 'hetero_tree_labels.txt') if l.strip()) != len(rmap): rmap = None
    for l in open(Path(T) / 'hetero_tree_leaves.txt'):
        w = l.split()
        if len(w) < 8 or w[0].startswith('#'): continue
        rows = np.loadtxt(w[7], dtype=int, ndmin=1)
        if rmap is not None: rows = rmap[rows - 1]
        tag = '%s__%s' % (name, w[0])
        if os.path.exists(w[5]) and os.path.exists(w[6]): spec.append('%s %s %s' % (tag, w[5], w[6]))
        leaves.append(dict(tree=name, id=w[0], kind=w[1], n=len(rows), nref=int(ref[rows].sum()), tag=tag,
                           depth=len(w[0].split('_')[0]) - 1))
(OUT / 'spec.txt').write_text('\n'.join(spec) + '\n')
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
subprocess.run([PY, str(SC / 'make_view.py'), str(OUT / 'views'), str(OUT / 'spec.txt'), '--box', '160'], check=True, env=env,
               stdout=subprocess.DEVNULL)
subprocess.run([PY, str(SC / 'arm_readout.py'), str(R160 / 'toy_consensus.mrc'), str(CL), str(OUT / 'arm.txt')] +
               [str(OUT / 'views' / (s.split()[0] + '.mrc')) for s in spec], check=True, env=env, stdout=subprocess.DEVNULL)
arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(OUT / 'arm.txt') if not l.startswith('#')}
for name, _ in TREES:
    lv = [x for x in leaves if x['tree'] == name]; tot = sum(x['nref'] for x in lv)
    for x in lv:
        x['rec'] = x['nref'] / max(tot, 1); x['pur'] = x['nref'] / max(x['n'], 1); x['f1'] = 2 * x['nref'] / max(x['n'] + tot, 1)
        x['ca'], x['oa'] = arm.get(x['tag'], (np.nan, np.nan))
    lv.sort(key=lambda x: -x['f1'])
    print('== %s: %d leaves, %d rows, reference closed rows in leaves %d (closed corner = closed-arm >= 0.95, open-arm <= 0.6)' %
          (name, len(lv), sum(x['n'] for x in lv), tot))
    b = lv[0]
    print('best closed-leaf F1 %.3f: %s (%d rows, ref %d, recall %.2f, purity %.3f, arm %.2f/%.2f)' %
          (b['f1'], b['id'], b['n'], b['nref'], b['rec'], b['pur'], b['ca'], b['oa']))
    cc = [x for x in lv if x['ca'] >= 0.95 and x['oa'] <= 0.6]
    print('closed-corner leaves: %d%s' % (len(cc), ''.join('\n    %s %d rows, ref %d (recall %.2f, purity %.3f, F1 %.3f), arm %.2f/%.2f' %
          (x['id'], x['n'], x['nref'], x['rec'], x['pur'], x['f1'], x['ca'], x['oa']) for x in cc)))
    up = sorted(lv, key=lambda x: -x['pur']); n = r = 0; best = (0., 0, 0, 0)
    for k, x in enumerate(up, 1):
        n += x['n']; r += x['nref']; f = 2 * r / max(n + tot, 1)
        if f > best[0]: best = (f, k, n, r)
    print('best union of purity-ranked leaves: F1 %.3f over %d leaves (%d rows, recall %.2f, purity %.3f): %s' %
          (best[0], best[1], best[2], best[3] / max(tot, 1), best[3] / max(best[2], 1), ' '.join(x['id'] for x in up[:best[1]])))
    print('%-14s %3s %-5s %7s %5s %6s %7s %6s %6s %6s' % ('leaf', 'd', 'kind', 'rows', 'ref', 'recall', 'purity', 'F1', 'closed', 'open'))
    for x in lv:
        print('%-14s %3d %-5s %7d %5d %6.3f %7.4f %6.3f %6.2f %6.2f' % (x['id'], x['depth'], x['kind'], x['n'], x['nref'], x['rec'],
              x['pur'], x['f1'], x['ca'], x['oa']))
