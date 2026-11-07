#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Which alignment are the root site candidates' maps (replay motion/cK/m.simple) in, and is each site's low-corr bin (1 vs 5) the same
# under the stage-B corr, m33's node corr (n.simple) and m.simple's own corr?
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import numpy as np
Q = '/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/'; G = '/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/gate/m33_root_m79/1_hetero_tree/node_r'
r = np.loadtxt(G + '/rows.txt', dtype=int) - 1
D = {'replay m.simple': np.load(Q + 'm79_root_old_1792953007/dump.npz'), 'replay n.simple': np.load(Q + 'm79_root/dump.npz'), 'replay b.simple': np.load(Q + 'm79_root_rows/dump.npz')}
P = {'m33 n.simple': np.load(Q + 'm79_root_rows/dumpP.npz'), 'm33 b.simple': np.load(Q + 'm33_root/dump.npz')}
for k in P: D[k] = {f: P[k][f][r] for f in ('state', 'corr', 'dfx', 'e1', 'e2')}
M = D['replay m.simple']
for k in D:
    if k == 'replay m.simple': continue
    X = D[k]; c = np.isfinite(M['corr']) & np.isfinite(X['corr'])
    print('m.simple vs %-15s pose median |de1| %5.1f |de2| %5.1f deg; corr r %.3f; dfx max |diff| %.1e' % (k, np.nanmedian(np.abs(((M['e1'] - X['e1']) + 180) % 360 - 180)),
          np.nanmedian(np.abs(M['e2'] - X['e2'])), np.corrcoef(M['corr'][c], X['corr'][c])[0, 1], np.nanmax(np.abs(M['dfx'] - X['dfx']))))
print('site  median corr bin1 / bin5 (lower side)   stage-B (replay b) | m33 n.simple | replay m.simple')
for s in (1, 2, 3, 5, 6, 7):
    b = np.loadtxt(G + '/motion/c%d/bins.txt' % s)
    out = []
    for k in ('replay b.simple', 'm33 n.simple', 'replay m.simple'):
        c = D[k]['corr']; m1, m5 = np.nanmedian(c[b == 1]), np.nanmedian(c[b == 5]); out.append('%.4f / %.4f (%s)' % (m1, m5, 'bin1' if m1 < m5 else 'bin5'))
    print('site %d  %s' % (s, ' | '.join(out)))
PYEOF
