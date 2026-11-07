#!/bin/bash
#SBATCH -J land_occ
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_occ.%j.out
# occupancy landscape in physical units (closed-site occupancy) from held-out windows + per-particle deconvolution. env: KIND e1b|real, RUN, [DUMP]
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
KIND=${KIND:?}; RUN=${RUN:?}; OUT=${OUT:-$RUN/occ}; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
if [ $KIND = e1b ]; then D="--sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 --poses $K_/lsamp/s30a/poses_combined.txt --truthmaps $K_/landscape/e0_s30/truth/leaves.txt --ctfabs 0 --mskdiam 230"
else D="--sim ${H}_A,${H}_B --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --cache $K_/lsamp/caches/rnasep692 --poses ${DUMP:?} --ctfabs 1 --mskdiam 232"; fi
$PY -u $T/land_occ.py --run $RUN --kind $KIND $D --cons $T/maps160/toy_consensus.mrc --closed $CL --out $OUT ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_occ $RUN DONE"
