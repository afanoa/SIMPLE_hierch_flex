#!/bin/bash
#SBATCH -J klem_lsamp
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape sampling, step 2 (lsamp.py) on a score space made by lsamp_score.sbatch. env SC (score dir), OUT, SIM, MR, CP, POSES, CTFABS,
# ARGS (further lsamp.py flags, e.g. "--naxes 1 --use r22_ --gt $T/simswing30/gt_particles.txt --codes swing:200:230"; commas allowed: set in
# the submitting shell)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
mkdir -p ${OUT:?}; echo "[$(date '+%F %T')] lsamp.py $ARGS -> $OUT"
$PY -u $T/lsamp.py --scores ${SC:?} --sim ${SIM:?} --mr ${MR:?} --cache ${CP:?} --poses ${POSES:?} --ctfabs ${CTFABS:-0} --out $OUT $ARGS 2>&1 | grep --line-buffered -v "Warning\|warn(" | tee $OUT/lsamp.log
echo "[$(date '+%F %T')] exit ${PIPESTATUS[0]}"
