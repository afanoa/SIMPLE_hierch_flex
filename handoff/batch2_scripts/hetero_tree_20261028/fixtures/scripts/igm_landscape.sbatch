#!/bin/bash
#SBATCH -J igm_land
#SBATCH -p norm
#SBATCH --constraint=x6342
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
# IgM landscapes (ole-repro's landscape_score.py) in the fulltoy_m69M root embedding: the given trees and the truth (genuine
# particles labelled combo + 1, junk 0), cut at depth 4 and at the leaves; env OUT, then NAME=TREE_DIR args
set -euo pipefail
: ${OUT:?}
PY=/mnt/nasapps/production/python/3.12/bin/python3; H=/mnt/beegfs/elmlund/afan/hetero_tree_20261028
E=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/fulltoy_m69M_igm/1_hetero_tree/node_r/flex/flex_pca_coordinates.txt
TR=$(dirname $OUT)/igm_truth_flat.txt
test -e $TR || $PY - $TR <<'PYEOF'
import sys, numpy as np
t = np.loadtxt('/mnt/beegfs/elmlund/afan/rnasep_klem/factorial/armsigm1_truth/arms_truth.txt', comments='#', ndmin=2)
lab = np.zeros(int(t[:, 0].max()), int); g = t[:, 2] == 0; lab[t[g, 0].astype(int) - 1] = t[g, 1].astype(int) + 1
np.savetxt(sys.argv[1], lab, fmt='%d')
PYEOF
export PYTHONDONTWRITEBYTECODE=1
{ echo "== depth 4"; $PY $H/ole_r14/tools/landscape_score.py --emb $E --depth 4 "$@" truth=$TR
  echo "== leaves"; $PY $H/ole_r14/tools/landscape_score.py --emb $E --leaves "$@" truth=$TR; } > $OUT 2>&1
