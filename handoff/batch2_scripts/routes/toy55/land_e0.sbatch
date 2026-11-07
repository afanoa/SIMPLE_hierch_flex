#!/bin/bash
#SBATCH -J land_e0
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_e0.%j.out
# oracle ceiling: swing30 particles scored against the 31 true swing maps (true and refined poses), then known-landscape recovery
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=$K_/landscape/e0_s30
mkdir -p $O; log(){ echo "[$(date '+%F %T') land_e0] $*"; }
[ -s $O/truth/labels.txt ] || $PY $T/land_truthmaps.py $T/maps160/swing30 $T/simswing30/gt_particles.txt $O/truth --ids 200:230 --mskdiam 200 || exit 1
for v in true node; do
  P=$T/simswing30/all.txt; [ $v = node ] && P=$K_/lsamp/s30a/poses_combined.txt
  [ -s $O/sc_$v/scores.npy ] || $PY -u $T/lsamp_score2.py --sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 --poses $P \
      --leaves $O/truth/leaves.txt --labels $O/truth/labels.txt --out $O/sc_$v --em_lim ${EMLIM:-10} 2>&1 | grep -v "Warning\|warn("
  log "scored ($v poses)"
  $PY -u $T/land_e0.py --scores $O/sc_$v --gt $T/simswing30/gt_particles.txt --out $O/ana_$v 2>&1 | grep -v "Warning\|warn("
done
log DONE
