#!/usr/bin/env python3
"""where does each tree first separate the RNase P closed state? Per finished node: the reference closed rows (8,720, the old
full692 tree's closed leaf) in each split group (count, recall of the node's, purity), the groups' stage-B state composition
(tree.log), and the closed-arm / open-arm readout of every stage-B state map (make_view 10 A views + arm_readout). Flags a node
whose split separates closed: a group holding >= 50 % of the node's reference rows at >= 2x the node's share, or stage-B state
maps whose closed-arm spread is >= 0.2 with the top state >= 0.8.
A row-local fixture run (its folder holds fixture_rows.txt and its root rows index the fixture) is mapped to the project rows; CN_NODES=id[,id]
limits the readout to those nodes. usage: closed_nodes2.py OUT_DIR NAME=TREE_DIR [NAME=TREE_DIR ..]"""
import sys, os, re, subprocess, numpy as np
from pathlib import Path
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
TREES = [a.split('=', 1) for a in sys.argv[2:]]
SC = Path(__file__).resolve().parent; PY = sys.executable; NR = 932496
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160')
CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
ref = np.zeros(NR + 1, bool)
ref[np.loadtxt('/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt', dtype=int)] = True
spec, nodes = [], []
ONLY = set(os.environ.get('CN_NODES', '').split(',')) - {''}
for name, T in TREES:
    T = Path(T); comp = {}; rmap, nrt = None, NR
    fx = T.parent / 'fixture_rows.txt'
    if fx.exists():
        m = np.loadtxt(fx, dtype=int, ndmin=1)
        if np.loadtxt(T / 'node_r' / 'rows.txt', dtype=int, ndmin=1).max() <= len(m): rmap, nrt = m, len(m)
    for l in open(T.parent / 'tree.log', errors='replace'):
        m = re.search(r'node (\S+): split group (\d+) \((\d+) particles; stage-B states ([0-9: ]+)\)', l)
        if m: comp[(m.group(1), int(m.group(2)))] = m.group(4).strip()
    for l in open(T / 'hetero_tree_nodes.txt'):
        w = l.split()
        if len(w) < 3 or int(w[2]) != 13: continue
        nid, depth = w[0], int(w[1]); d = T / ('node_' + nid)
        if ONLY and nid not in ONLY: continue
        rows = np.loadtxt(d / 'rows.txt', dtype=int, ndmin=1)
        f = d / 'split/split_children.txt'
        if not f.exists(): continue
        ch = np.loadtxt(f, dtype=int, ndmin=1)
        if len(ch) == nrt: ch = ch[rows - 1]
        if len(ch) != len(rows): continue
        r = ref[rows if rmap is None else rmap[rows - 1]]; groups = []
        for g in range(1, ch.max() + 1):
            m = ch == g
            if m.sum(): groups.append((g, int(m.sum()), int((m & r).sum()), comp.get((nid, g), '')))
        states = []
        for e in sorted((d / 'stageB').glob('recvol_state??_even_unfil.mrc')):
            k = e.name[12:14]; o = e.with_name(e.name.replace('_even_', '_odd_'))
            if o.exists():
                tag = '%s__%s__s%s' % (name, nid, k); spec.append('%s %s %s' % (tag, e, o)); states.append((k, tag))
        nodes.append(dict(tree=name, id=nid, depth=depth, n=len(rows), nref=int(r.sum()), groups=groups, states=states))
(OUT / 'spec.txt').write_text('\n'.join(spec) + '\n')
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
subprocess.run([PY, str(SC / 'make_view.py'), str(OUT / 'views'), str(OUT / 'spec.txt'), '--box', '160'], check=True, env=env,
               stdout=subprocess.DEVNULL)
subprocess.run([PY, str(SC / 'arm_readout.py'), str(R160 / 'toy_consensus.mrc'), str(CL), str(OUT / 'arm.txt')] +
               [str(OUT / 'views' / (s.split()[0] + '.mrc')) for s in spec], check=True, env=env, stdout=subprocess.DEVNULL)
arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(OUT / 'arm.txt') if not l.startswith('#')}
for name, _ in TREES:
    print('== %s (reference closed rows: 8,720; readout closed-arm / open-arm; closed corner >= 0.95 / <= 0.6)' % name)
    for nd in sorted([x for x in nodes if x['tree'] == name], key=lambda x: (x['depth'], x['id'])):
        sh = nd['nref'] / max(nd['n'], 1)
        flag_ref = any(g[2] >= 0.5 * nd['nref'] and g[2] / max(g[1], 1) >= 2 * sh for g in nd['groups']) and nd['nref'] >= 500
        ca = [arm[t][0] for _, t in nd['states'] if t in arm]
        flag_map = len(ca) >= 2 and max(ca) - min(ca) >= 0.2 and max(ca) >= 0.8
        print('%-9s d%d %7d rows | ref %5d (share %.4f) | %s%s' % (nd['id'], nd['depth'], nd['n'], nd['nref'], sh,
              'SEPARATES(ref) ' if flag_ref else '', 'SEPARATES(map)' if flag_map else ''))
        for g, n, nr, c in nd['groups']:
            print('    group %d %7d rows | ref %5d (recall %.2f, purity %.4f) | stage-B states %s' % (g, n, nr, nr / max(nd['nref'], 1),
                  nr / max(n, 1), c))
        print('    stage-B maps: ' + '  '.join('s%s %.2f/%.2f' % (k, *arm[t]) for k, t in nd['states'] if t in arm))
