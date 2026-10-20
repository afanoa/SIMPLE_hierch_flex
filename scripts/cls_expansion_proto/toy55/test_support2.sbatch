#!/bin/bash
#SBATCH -J klem_support_test2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# (a) truth readout of the toy state sets from the support test (off vs on); (b) the real boards with the cap lifted
# (maxk 30), support off vs on: how many states the data support per board.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; O=$A/rnasep_klem/minpop_test; GT=$W/simah/gt_particles.txt
echo "######## (a) toy truth readout"
for b in toyah_tree_fam1_child toyah_tree_fam2_child; do for m in off on; do $PY $W/score_states.py $O/${b}_$m/states_clusters.txt $GT x "$b $m"; done; done
echo "######## (b) real boards, cap lifted"
for b in f692/fam1/region_1/board f692/fam1/regionmotion_1/board f692/fam2/region_1/board f692/fam2/regionmotion_1/board f692/fam1/klem_fb; do tag=$(echo $b | tr / _)
  for m in off on; do echo "=================== $b  support $m  maxk 30"; $PY $W/node_states2.py $A/rnasep_klem/$b $O/${tag}_k30_$m --nrows 932496 --sim 0.7 --nmax 30000 --maxk 30 $([ $m = on ] && echo --support) 2>&1 | grep -v Warning | grep -v "^  state"; done; done
echo "ALL DONE"
