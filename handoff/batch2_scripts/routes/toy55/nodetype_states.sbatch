#!/bin/bash
#SBATCH -J nodetype_states
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Hard-state geometry (nodetype_states.py) of every node's focused refinement. env SET: s30a (truth) | real (f692rf, no truth)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
case ${SET:-s30a} in
  s30a) R=$K_/rflex_s30a; for u in r r1 r2 r11 r12 r21 r22; do $PY -u $T/nodetype_states.py $(ls -td $R/node_$u/*_refine3D | head -1) $R/node_$u/refined_states.txt --gt $T/simswing30/gt_particles.txt --tmaps $T/simswing30/maps --tag "[s30a $u]" 2>&1 | grep --line-buffered -v Warning; done ;;
  real) R=$K_/rflex_f692rf; for u in r r1 r2 r21 r22; do $PY -u $T/nodetype_states.py $(ls -td $R/node_$u/*_refine3D | head -1) $R/node_$u/refined_states.txt --mskdiam 232 --tag "[f692rf $u]" 2>&1 | grep --line-buffered -v Warning; done ;;
esac
echo "[$(date '+%F %T') nodetype_states ${SET:-s30a}] DONE"
