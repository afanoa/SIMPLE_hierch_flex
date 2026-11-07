#!/bin/bash
#SBATCH -J subset_only
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# no-op CHAIN for synth_tree_prep_full: writes only the parent subset (every flagged row, 1-based) of the run dir, which
# dtree_v10.sbatch copies to its root node (env RTAG, METAALL from synth_tree_prep_full)
PY=/mnt/nasapps/production/python/3.12/bin/python3; R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}
$PY -c "import numpy as np; f = np.load('${METAALL:?}')['flags'].astype(int); rows = np.where(f == 1)[0] + 1; np.savetxt('$R/subset_all.txt', rows, fmt='%d'); print('$RTAG parent subset: %d of %d rows' % (len(rows), len(f)))"
