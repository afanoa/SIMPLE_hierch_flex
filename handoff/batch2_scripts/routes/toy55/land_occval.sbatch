#!/bin/bash
#SBATCH -J land_occval
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_occval.%j.out
# synthetic validation of the occupancy route: land_lam (coordinate) then land_occ (closed-site occupancy units) on one imposed shape. env SET SH
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
SET=${SET:-fam}; SH=${SH:?}; O=$K_/landscape/occval_s30_$SET/$SH; mkdir -p $O; USE=r22_; FAM=0; [ $SET = fam ] && { USE=; FAM=1; }
[ -s $O/e1b_particles.npz ] || $PY -u $T/land_lam.py --scores $K_/lsamp/s30a_ps2 --use "$USE" --family $FAM --sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 \
   --poses $K_/lsamp/s30a/poses_combined.txt --gt $T/simswing30/gt_particles.txt --codes swing:200:230 --shape $SH --out $O 2>&1 | grep --line-buffered -v "Warning\|warn("
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY -u $T/land_occ.py --run $O --kind e1b --sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 --poses $K_/lsamp/s30a/poses_combined.txt \
   --truthmaps $K_/landscape/e0_s30/truth/leaves.txt --ctfabs 0 --mskdiam 230 --cons $T/maps160/toy_consensus.mrc --closed $CL --out $O/occ 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_occval $SET $SH DONE"
