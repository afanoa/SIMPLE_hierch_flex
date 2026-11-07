#!/usr/bin/env python3
"""one node's stage-B states and delivered split groups side by side across runs: per stage-B state (stageB/prob_state_evidence.txt
line j = stageB/node_rows.txt row j) its size, reference closed rows (the old full692 closed leaf, 8,720), purity and the closed-arm /
open-arm readout of its map (make_view 10 A, box 160 + arm_readout); per delivered group (split/split_children.txt) its size, make-up
by stage-B state, reference rows, purity and F1 against the node's and all 8,720 reference rows; K*. A row-local fixture run (its
folder holds fixture_rows.txt and its root rows index the fixture) is mapped to the project rows.
usage: compare_stageb.py OUT_DIR NODE NAME=TREE_DIR [NAME=TREE_DIR ..]"""
import sys, os, subprocess, numpy as np
from pathlib import Path
OUT, NODE = Path(sys.argv[1]), sys.argv[2]; OUT.mkdir(parents=True, exist_ok=True)
SC = Path(__file__).resolve().parent; PY = sys.executable; NR = 932496
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160'); CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
ref = np.zeros(NR + 1, bool)
ref[np.loadtxt('/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt', dtype=int)] = True
TOT = int(ref.sum()); runs, spec = [], []
for arg in sys.argv[3:]:
    name, T = arg.split('=', 1); T = Path(T); d = T / ('node_' + NODE)
    rmap = None; fx = T.parent / 'fixture_rows.txt'
    if fx.exists():
        m = np.loadtxt(fx, dtype=int, ndmin=1)
        if np.loadtxt(T / 'node_r' / 'rows.txt', dtype=int, ndmin=1).max() <= len(m): rmap = m
    glob = lambda r: r if rmap is None else rmap[r - 1]
    nr = np.loadtxt(d / 'stageB' / 'node_rows.txt', dtype=int, ndmin=1)[1:]
    ev = [l.split() for l in open(d / 'stageB' / 'prob_state_evidence.txt') if l.strip() and not l.startswith('#')]
    st = np.array([int(w[1]) for w in ev]); assert len(st) == len(nr), (name, len(st), len(nr))
    gr = glob(nr); rr = ref[gr]
    rows = np.loadtxt(d / 'rows.txt', dtype=int, ndmin=1); ch = np.loadtxt(d / 'split' / 'split_children.txt', dtype=int, ndmin=1)
    nproj = len(rmap) if rmap is not None else NR
    chn = ch[nr - 1] if len(ch) == nproj else (ch if len(ch) == len(nr) else None)
    states = {}
    for s in sorted(np.unique(st)):
        if s <= 0: continue
        e = d / 'stageB' / ('recvol_state%02d_even_unfil.mrc' % s); o = d / 'stageB' / ('recvol_state%02d_odd_unfil.mrc' % s)
        tag = '%s__s%02d' % (name, s)
        if e.exists() and o.exists(): spec.append('%s %s %s' % (tag, e, o))
        states[s] = dict(n=int((st == s).sum()), nref=int(rr[st == s].sum()), tag=tag)
    groups = []
    if chn is not None:
        for g in sorted(np.unique(chn)):
            m = chn == g
            groups.append(dict(g=int(g), n=int(m.sum()), nref=int(rr[m].sum()),
                               comp=' '.join('%d:%d' % (s, int(((st == s) & m).sum())) for s in sorted(states) if ((st == s) & m).any())))
    runs.append(dict(name=name, n=len(nr), nref=int(rr.sum()), states=states, groups=groups, k=len([g for g in groups if g['g'] > 0])))
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1'); arm = {}
if spec:
    (OUT / 'spec.txt').write_text('\n'.join(spec) + '\n')
    subprocess.run([PY, str(SC / 'make_view.py'), str(OUT / 'views'), str(OUT / 'spec.txt'), '--box', '160'], check=True, env=env, stdout=subprocess.DEVNULL)
    subprocess.run([PY, str(SC / 'arm_readout.py'), str(R160 / 'toy_consensus.mrc'), str(CL), str(OUT / 'arm.txt')] +
                   [str(OUT / 'views' / (s.split()[0] + '.mrc')) for s in spec], check=True, env=env, stdout=subprocess.DEVNULL)
    arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(OUT / 'arm.txt') if not l.startswith('#')}
print('node %s (reference closed rows: %d in all; closed corner = closed-arm >= 0.95, open-arm <= 0.6)' % (NODE, TOT))
for r in runs:
    print('== %s: %d rows, %d reference rows (share %.4f), K* %d' % (r['name'], r['n'], r['nref'], r['nref'] / max(r['n'], 1), r['k']))
    for s, x in r['states'].items():
        a = arm.get(x['tag'])
        print('   stage-B state %d: %6d rows | ref %5d, purity %.4f | map %s' % (s, x['n'], x['nref'], x['nref'] / max(x['n'], 1),
              '%.2f / %.2f' % a if a else 'not found'))
    for g in r['groups']:
        lab = 'junk ' if g['g'] == 0 else 'group %d' % g['g']
        print('   %-8s %6d rows | stage-B make-up %s | ref %5d, purity %.4f, F1 node %.3f, F1 all %.3f' % (lab, g['n'], g['comp'], g['nref'],
              g['nref'] / max(g['n'], 1), 2 * g['nref'] / max(g['n'] + r['nref'], 1), 2 * g['nref'] / max(g['n'] + TOT, 1)))
