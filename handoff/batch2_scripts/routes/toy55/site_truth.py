#!/usr/bin/env python3
"""Per-row site truth for the factorial-site tests (master row order): simsites (combination id -> domain A present, blob B present, arm C
closed) or simarmsigm1 (conformer -> 5 arms x theta / phi; binary site states = the sign of each angle). Output columns:
row combo junk then one column per site variable; header '# row combo junk <names>' and a line '# kind <names: d (discrete) | c>'.
usage: site_truth.py sites GT COMBOS OUT | site_truth.py arms GT TRUTH_ARMS OUT"""
import sys, numpy as np
mode, gt, aux, out = sys.argv[1:5]
G = np.loadtxt(gt, comments='#'); N = int(G[:, 0].max()); conf = np.zeros(N, int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int); junk = (conf == 0).astype(int)
if mode == 'sites':
    C = {int(r[0]): r[1:4].astype(int) for r in np.loadtxt(aux, comments='#')}; names = ['A_domain', 'B_blob', 'C_closed']; kinds = 'ddd'
    V = np.array([C[k] if k > 0 else [-1, -1, -1] for k in conf], float)
else:
    T = {int(w[0]): [float(x) for x in w[1:]] for w in (l.split() for l in open(aux)) if w and not w[0].startswith('#')}; na = len(next(iter(T.values()))) // 2
    names = ['a%d_%s' % (k + 1, t) for k in range(na) for t in ('th', 'ph')]; kinds = 'c' * len(names)
    V = np.array([T[k] if k > 0 else [np.nan] * len(names) for k in conf], float)
with open(out, 'w') as f:
    f.write('# row combo junk %s\n# kind %s\n' % (' '.join(names), ' '.join(kinds)))
    for i in range(N): f.write('%d %d %d %s\n' % (i + 1, conf[i], junk[i], ' '.join(('%d' % v if kinds[j] == 'd' else '%.2f' % v) if np.isfinite(v) else 'nan' for j, v in enumerate(V[i]))))
print('%s: %d rows, %d junk, sites %s' % (out, N, junk.sum(), ' '.join(names)))
