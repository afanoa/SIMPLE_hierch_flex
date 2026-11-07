#!/bin/bash
#SBATCH -J land_pscore
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_pscore.%j.out
# per-particle path likelihoods with local pose integration. env MAPS OUT KIND (s30|real) [ROWS POSES ARGS]
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
if [ ${KIND:-s30} = s30 ]; then D="--sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 --poses ${POSES:-$K_/lsamp/s30a/poses_combined.txt} --rowflags ${ROWS:-$K_/landscape/e0_s30/truth/labels.txt} --ctfabs 0 --mskdiam 230"
else D="--sim ${H}_A,${H}_B --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --cache $K_/lsamp/caches/rnasep692 --poses ${POSES:?} --rowflags ${ROWS:?} --ctfabs 1 --mskdiam 232"; fi
$PY -u $T/land_pscore.py $D --maps ${MAPS:?} --out ${OUT:?} ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_pscore DONE"
