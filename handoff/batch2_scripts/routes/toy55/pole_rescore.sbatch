#!/bin/bash
#SBATCH -J klem_pole_rescore
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# null_score3.py (no pole-rule index) on the 14 nodes whose seeding score the null and RNase P series used, compared with the score
# file each series binned on: the largest abs difference over the rows
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port; E=$K/factorial/phaseE
IG=$S/fulltoy_m52M_igm/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; M=$S/full692_m33_s1/1_hetero_tree; W=$E/pole_rescore; mkdir -p $W
for e in $CD/node_r111:$E/null/null_cdisc_r111/work/score.txt $CD/node_r1111:$E/null/null_cdisc_r1111/work/score.txt $IG/node_r212:$E/null/null_igm_r212/work/score.txt \
         $IG/node_r2122:$E/null/null_igm_r2122/work/score.txt $CD/node_r11111:$E/null2/null_cdisc_r11111/work/score.txt \
         $IG/node_r1111112:$E/null3/null_igm_r1111112/work/score.txt $IG/node_r111212:$E/null3/null_igm_r111212/work/score.txt \
         $IG/node_r12121:$E/null3/null_igm_r12121/work/score.txt $IG/node_r211112:$E/null3/null_igm_r211112/work/score.txt \
         $CD/node_r1112:$E/null3/null_cdisc_r1112/work/score.txt $CD/node_r111111:$E/null3/null_cdisc_r111111/work/score.txt \
         $M/node_r:$E/rnasep_m33/node_r/scores/pole.txt $M/node_r1:$E/rnasep_m33/node_r1/scores/pole.txt $M/node_r2:$E/rnasep_m33/node_r2/scores/pole.txt; do
  N=${e%%:*}; old=${e#*:}; nm=$(basename $N); new=$W/$nm.txt
  msg=$($PY $T/null_score3.py $N $new 2>&1 | tail -1)
  d=$($PY -c "import numpy as np, sys; a = np.loadtxt(sys.argv[1]); b = np.loadtxt(sys.argv[2]); print('%.2e' % np.abs(a - b).max() if a.shape == b.shape else 'shape %s vs %s' % (a.shape, b.shape))" $new $old 2>&1 | tail -1)
  echo "$nm | max abs diff vs the series' score $d | $msg"
done
