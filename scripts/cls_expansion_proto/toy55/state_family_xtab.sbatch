#!/bin/bash
#SBATCH -J klem_famxtab
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# cross-tab of merged states against family membership (env MODES, XS)
PY=/mnt/nasapps/production/python/3.12/bin/python3; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for M in ${MODES:-null_50k_A null_50k_B two_50k_A two_50k_B}; do R=$K/synth_$M; echo "== $M ($XS)"
$PY - <<PYEOF
import numpy as np
st = np.loadtxt('$R/consolidated_xtree$XS/states.txt').astype(int); fam = np.zeros(len(st), int)
for f in (1, 2):
    rows = np.loadtxt('$R/expand/family%d_subset.txt' % f).astype(int) - 1; fam[rows[rows < len(st)]] = f
for s in np.unique(st[st > 0]): m = st == s; print('  state %d: %6d particles, family 1 %.3f, family 2 %.3f' % (s, m.sum(), np.mean(fam[m] == 1), np.mean(fam[m] == 2)))
PYEOF
done
