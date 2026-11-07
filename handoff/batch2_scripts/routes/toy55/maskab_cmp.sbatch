#!/bin/bash
#SBATCH -J maskab_cmp
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Three-way comparison (original current-mask run, current-mask replicate, children-mask run) at r22 and r12: labels vs truth and
# pairwise agreement (maskab_compare.py), plus the per-state map summaries already printed by the maskab jobs.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
for u in r22 r12; do echo "===== node $u"
  $PY -u $T/maskab_compare.py $T/simswing30/gt_particles.txt original_current=$K_/rflex_s30a/node_$u/refined_states.txt replicate_current=$K_/maskab_s30a_${u}_mask/refined_states.txt children_mask=$K_/maskab_s30a_${u}_kids_mask/refined_states.txt 2>&1 | grep -v Warning
  for f in $(grep -l "maskab $u\]" $K_/maskab.*.out); do grep -h "member-weighted mean\|refine3D .* s," $f | sed "s#^#  $(basename $f): #" | cut -c1-260; done
done
