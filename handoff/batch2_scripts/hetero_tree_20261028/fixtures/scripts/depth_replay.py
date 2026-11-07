#!/usr/bin/env python3
"""depth cap replay from node rows only: the tree cut at depth D puts every particle in its ancestor at D (a lineage that ended
shallower keeps its last node). Per depth: groups, median size, best closed F1 (recall, purity) against the reference closed rows
(8,720, the old full692 closed leaf), the best union of purity-ranked groups, and the closed-arm / open-arm readout of the best
group's map (its parent's stage-B map of the group's dominant state; at the leaf row, the leaf map), where that map was kept.
usage: depth_replay.py OUT_DIR NAME=TREE_DIR:D1,D2,.. [..]   (TREE_DIR = 1_hetero_tree; leaf row added when hetero_tree_leaves.txt exists)"""
import sys, os, re, subprocess, numpy as np
from pathlib import Path
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
SC = Path(__file__).resolve().parent; PY = sys.executable; NR = 932496
SCR, BEE = '/scratch/cluster_scratch/atifao', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/from_scratch'
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160'); CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
ref = np.zeros(NR + 1, bool)
ref[np.loadtxt('/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt', dtype=int)] = True
def kept(p):   # the beegfs copy of a scratch path, else the scratch file while it exists
    p = str(p); b = p.replace(SCR, BEE)
    return b if os.path.exists(b) else (p if os.path.exists(p) else None)
def score(groups, tot):   # groups: [(gid, rows)] -> per-group stats, best F1, best purity-ranked union
    st = [(g, len(r), int(ref[r].sum())) for g, r in groups]
    best = max(st, key=lambda x: 2 * x[2] / (x[1] + tot))
    up = sorted(st, key=lambda x: -x[2] / max(x[1], 1)); n = c = 0; bu = (0., 0, 0, 0)
    for k, x in enumerate(up, 1):
        n += x[1]; c += x[2]; f = 2 * c / (n + tot)
        if f > bu[0]: bu = (f, k, n, c)
    return st, best, bu
rows_out, spec = [], []
for arg in sys.argv[2:]:
    name, rest = arg.split('=', 1); T, ds = rest.rsplit(':', 1); T = Path(T); depths = [int(x) for x in ds.split(',')]
    tlog = T.parent / 'tree.log'; comp = {}
    if tlog.exists():
        for l in open(tlog, errors='replace'):
            m = re.search(r'node (\S+): split group (\d+) \((\d+) particles; stage-B states ([0-9: ]+)\)', l)
            if m: comp[(m.group(1), int(m.group(2)))] = [tuple(map(int, s.split(':'))) for s in m.group(4).split()]
    ids = [l.split()[0] for l in open(T / 'hetero_tree_nodes.txt') if l.split() and not l.startswith('#')]
    have = {i for i in ids if (T / ('node_' + i) / 'rows.txt').exists()}
    rows = {i: np.loadtxt(T / ('node_' + i) / 'rows.txt', dtype=int, ndmin=1) for i in have}
    root = rows['r']; tot = int(ref[root].sum())
    for D in depths:
        g = [i for i in have if len(i) - 1 == D or (len(i) - 1 < D and not any(j != i and j.startswith(i) and len(j) == len(i) + 1 for j in have))]
        cov = sum(len(rows[i]) for i in g)
        st, best, bu = score([(i, rows[i]) for i in g], tot)
        gid = best[0]; p, k = gid[:-1], int(gid[-1]) if len(gid) > 1 else 0
        mp = None
        if p and (p, k) in comp:
            s = max(comp[(p, k)], key=lambda x: x[1])[0]
            e = kept(T / ('node_' + p) / 'stageB' / ('recvol_state%02d_even_unfil.mrc' % s)); o = e and kept(e.replace('_even_', '_odd_'))
            if e and o: mp = ('%s__d%d__%s_s%02d' % (name, D, p, s), e, o); spec.append('%s %s %s' % mp)
        rows_out.append(dict(tree=name, depth='%d' % D, ng=len(g), med=int(np.median([x[1] for x in st])), cov=cov, n=len(root),
                             best=best, tot=tot, bu=bu, map=mp[0] if mp else None, mapsrc='%s state %d' % (p, s) if mp else 'not kept'))
    lf = T / 'hetero_tree_leaves.txt'
    if lf.exists() and (T / 'hetero_tree_labels.txt').exists():
        lv = []
        for l in open(lf):
            w = l.split()
            if len(w) < 8 or w[0].startswith('#'): continue
            r = kept(w[7]); lv.append((w[0], np.loadtxt(r, dtype=int, ndmin=1), w[5], w[6]))
        st, best, bu = score([(x[0], x[1]) for x in lv], tot)
        b = [x for x in lv if x[0] == best[0]][0]; e, o = kept(b[2]), kept(b[3]); mp = None
        if e and o: mp = ('%s__leaves__%s' % (name, b[0]), e, o); spec.append('%s %s %s' % mp)
        rows_out.append(dict(tree=name, depth='leaves', ng=len(lv), med=int(np.median([x[1] for x in st])), cov=sum(x[1] for x in st),
                             n=len(root), best=best, tot=tot, bu=bu, map=mp[0] if mp else None, mapsrc='leaf map' if mp else 'not kept'))
arm = {}
if spec:
    (OUT / 'spec.txt').write_text('\n'.join(spec) + '\n'); env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([PY, str(SC / 'make_view.py'), str(OUT / 'views'), str(OUT / 'spec.txt'), '--box', '160'], check=True, env=env, stdout=subprocess.DEVNULL)
    subprocess.run([PY, str(SC / 'arm_readout.py'), str(R160 / 'toy_consensus.mrc'), str(CL), str(OUT / 'arm.txt')] +
                   [str(OUT / 'views' / (s.split()[0] + '.mrc')) for s in spec], check=True, env=env, stdout=subprocess.DEVNULL)
    arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(OUT / 'arm.txt') if not l.startswith('#')}
print('reference closed rows per tree in its root: ' + ', '.join(sorted({'%s %d' % (r['tree'], r['tot']) for r in rows_out})))
print('%-6s %-6s %6s %7s %9s | %-9s %7s %5s %6s %7s %6s | %-24s | %s' % ('tree', 'depth', 'groups', 'median', 'covered', 'best', 'rows',
      'ref', 'recall', 'purity', 'F1', 'best union F1 (k, recall, purity)', 'best-group map closed/open (source)'))
for r in rows_out:
    g, n, c = r['best']; f, k, un, uc = r['bu']; a = arm.get(r['map'])
    print('%-6s %-6s %6d %7d %9s | %-9s %7d %5d %6.3f %7.4f %6.3f | %.3f (%d, %.2f, %.3f)%s | %s' % (r['tree'], r['depth'], r['ng'], r['med'],
          '%d/%d' % (r['cov'], r['n']), g, n, c, c / r['tot'], c / n, 2 * c / (n + r['tot']), f, k, uc / r['tot'], uc / max(un, 1),
          ' ' * 2, ('%.2f / %.2f (%s)' % (a[0], a[1], r['mapsrc'])) if a else r['mapsrc']))
