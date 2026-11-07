#!/bin/bash
#SBATCH -J r11child
#SBATCH -p csbdevel
#SBATCH --constraint=x6342
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
# next level of the r11 test: one r11 arm's stage-B state (PARENT run dir, STATE) as a root-sized fixture FIX, then one sst_gate4 node in this job
set -euo pipefail
: ${PARENT:?} ${STATE:?} ${FIX:?}
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/closed_r11_routeA
S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
LR=$PARENT/1_hetero_tree/node_r/leaf_r_s${STATE}_rows.txt
if [ ! -e $FIX/manifest.txt ]; then
  test ! -e $FIX; mkdir -p $FIX; cp $P/input.simple $FIX/input.simple; cp $LR $FIX/fixture_rows.txt
  n=$(grep -c . $P/fixture_selection.txt)
  awk -v n=$n '{s[$1] = 1} END {for (r = 1; r <= n; r++) print ((r in s) ? 1 : 0)}' $LR > $FIX/fixture_selection.txt
  printf '%s\n' "source_node=$LR" "mode=root-sized" "rows=$(grep -c . $FIX/fixture_rows.txt)" "selected=$(awk '$1 > 0' $FIX/fixture_selection.txt | wc -l)" \
    "nrows=$n" "source_project_sha256=$(sha256sum $P/input.simple | cut -c1-64)" "made=$(date -Is)" > $FIX/manifest.txt
fi
export FIX; bash $S/sst_gate4.sbatch
