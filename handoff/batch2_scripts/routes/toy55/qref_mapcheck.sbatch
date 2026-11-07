#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Does the replay's node-rows-only stage-B dump map onto the full project through rows.txt (dfx per row must agree), and whose
# orientations / corr does it share: m33's stage B (b.simple) or m33's node refinement (n.simple)?
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import numpy as np
Q = '/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/'
A = np.load(Q + 'm79_root_rows/dump.npz'); P = np.load(Q + 'm79_root_rows/dumpP.npz'); M = np.load(Q + 'm33_root/dump.npz')
r = np.loadtxt('/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/gate/m33_root_m79/1_hetero_tree/node_r/rows.txt', dtype=int) - 1
print('rows %d, ascending %s; replay b %d, m33 n %d, m33 b %d' % (len(r), bool((np.diff(r) > 0).all()), len(A['dfx']), len(P['dfx']), len(M['dfx'])))
for nm, X in (('m33 n.simple', P), ('m33 b.simple', M)):
    d = np.abs(A['dfx'] - X['dfx'][r]); ok = np.isfinite(d)
    c = np.isfinite(A['corr']) & np.isfinite(X['corr'][r])
    e1 = np.abs(((A['e1'] - X['e1'][r]) + 180) % 360 - 180); e2 = np.abs(A['e2'] - X['e2'][r])
    print('%-13s dfx max |diff| %.2e (%d rows); corr r %.3f, medians %.4f vs %.4f; median |de1| %.1f, |de2| %.1f deg'
          % (nm, d[ok].max(), ok.sum(), np.corrcoef(A['corr'][c], X['corr'][r][c])[0, 1], np.nanmedian(A['corr']), np.nanmedian(X['corr'][r]), np.nanmedian(e1), np.nanmedian(e2)))
c = np.isfinite(P['corr'][r]) & np.isfinite(M['corr'][r])
print('m33 n vs m33 b: corr r %.3f; median |de1| %.1f, |de2| %.1f deg' % (np.corrcoef(P['corr'][r][c], M['corr'][r][c])[0, 1],
      np.nanmedian(np.abs(((P['e1'][r] - M['e1'][r]) + 180) % 360 - 180)), np.nanmedian(np.abs(P['e2'][r] - M['e2'][r]))))
PYEOF
