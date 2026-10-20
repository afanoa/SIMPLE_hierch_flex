#!/bin/bash
#SBATCH -J klem_toy_readout
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
echo "ALL DONE"
