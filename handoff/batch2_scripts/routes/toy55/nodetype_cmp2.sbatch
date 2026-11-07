#!/bin/bash
#SBATCH -J nodetype_cmp2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Hard states (assigned at lp 10) vs score-loop coordinate leaves (assigned at lp 8) at the s30a nodes that have both, in the band
# neither used (8-6 A), with truth (nodetype_compare.py).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_s30a; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for u in r11 r12 r21 r22; do
  $PY -u $T/nodetype_compare.py --a hard:$(ls -td $R/node_$u/*_refine3D | head -1):$R/node_$u/refined_states.txt --b coordleaves:$R/coord_$u/leaves/1_reconstruct3D:$R/coord_$u/leaf_states.txt \
    --band 8,6 --gt $T/simswing30/gt_particles.txt --tmaps $T/simswing30/maps --tag "[$u]" 2>&1 | grep --line-buffered -v Warning; done
echo "[$(date '+%F %T') nodetype_cmp2] DONE"
