#!/bin/bash
#SBATCH -J nodetype_cmp
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Hard states vs cheap coordinate bins per node (nodetype_compare.py) in bands neither assignment used. env SET s30a | real
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
case ${SET:-s30a} in
  s30a) R=$K_/rflex_s30a; NODES="r r1 r2 r11 r12 r21 r22"; MD=220; TRU="--gt $T/simswing30/gt_particles.txt --tmaps $T/simswing30/maps"; KD=cheapkids_s30a ;;
  real) R=$K_/rflex_f692rf; NODES="r r1 r2 r21 r22"; MD=232; TRU=""; KD=cheapkids_f692rf ;;
esac
for u in $NODES; do for band in 10,8 15,10; do
  $PY -u $T/nodetype_compare.py --a hard:$(ls -td $R/node_$u/*_refine3D | head -1):$R/node_$u/refined_states.txt --b bins:$(ls -td $K_/${KD}_$u/rec/*_reconstruct3D | head -1):$K_/${KD}_$u/children.txt \
    --coords $K_/${KD}_$u/coords.txt --band $band --mskdiam $MD $TRU --tag "[$u]" 2>&1 | grep --line-buffered -v Warning; done; done
echo "[$(date '+%F %T') nodetype_compare ${SET:-s30a}] DONE"
