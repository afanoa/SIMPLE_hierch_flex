#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# qtoy_score.py on two fake trees from the truth: "perfect" (one leaf per state, junk named) and "quality" (leaves by B tercile x state):
# recall / precision must read 1 / 1 for perfect, the quality eta^2 about 0 for perfect and high for quality
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; TR=$T/simqtoy/truth_qtoy.txt; W=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/qtoy_score_test
mkdir -p $W/perfect/1_hetero_tree $W/quality/1_hetero_tree; cd $W
$PY - $TR <<'PYEOF'
import sys, numpy as np
T = np.genfromtxt(sys.argv[1], comments='#'); conf, jt, B = T[:, 7].astype(int), T[:, 6].astype(int), T[:, 11]; st = np.where(jt > 0, -jt, conf)
u = sorted(set(st.tolist())); lab = np.array([u.index(s) + 1 for s in st]); np.savetxt('perfect/1_hetero_tree/hetero_tree_labels.txt', lab, fmt='%d')
open('perfect/1_hetero_tree/hetero_tree_leaves.txt', 'w').write(''.join('L%d %s %d\n' % (k + 1, 'junk' if s < 0 else 'leaf', (st == s).sum()) for k, s in enumerate(u)))
t = np.digitize(B, [66.7, 133.3]); lab2 = (lab - 1) * 3 + t + 1; np.savetxt('quality/1_hetero_tree/hetero_tree_labels.txt', lab2, fmt='%d')
open('quality/1_hetero_tree/hetero_tree_leaves.txt', 'w').write(''.join('L%d_%d %s %d\n' % (k + 1, j, 'junk' if s < 0 else 'leaf', ((st == s) & (t == j)).sum()) for k, s in enumerate(u) for j in range(3)))
PYEOF
for d in perfect quality; do echo "=== $d"; $PY $T/qtoy_score.py $W/$d/1_hetero_tree $TR | head -9 | cut -c1-600; done
