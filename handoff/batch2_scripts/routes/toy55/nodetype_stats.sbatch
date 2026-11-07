#!/bin/bash
#SBATCH -J nodetype_stats
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Node-type statistics (nodetype_stats.py) for a list of cheap-children dirs. env SET: s30a (truth, swing30 maps), dev (null/two r and
# rc, truth from maps160/docked), real (f692rf, no truth)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
case ${SET:-s30a} in
  s30a) for u in r r1 r2 r11 r12 r21 r22; do $PY -u $T/nodetype_stats.py $K_/cheapkids_s30a_$u --gt $T/simswing30/gt_particles.txt --tmaps $T/simswing30/maps --tag "[s30a $u]" 2>&1 | grep --line-buffered -v Warning; done ;;
  dev)  for m in null_50k two_50k; do for u in r rc; do $PY -u $T/nodetype_stats.py $K_/cheapkids_root_${m}_$u --gt $T/sim$m/gt_particles.txt --tmaps $T/maps160/docked --tag "[$m $u]" 2>&1 | grep --line-buffered -v Warning; done; done ;;
  real) for u in r r1 r2 r21 r22; do $PY -u $T/nodetype_stats.py $K_/cheapkids_f692rf_$u --mskdiam 232 --tag "[f692rf $u]" 2>&1 | grep --line-buffered -v Warning; done ;;
esac
echo "[$(date '+%F %T') nodetype_stats ${SET:-s30a}] DONE"
