#!/bin/bash
#SBATCH -J land_e1b
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 05:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_e1b.%j.out
# truth-free landscape tests on swing30 (s30a tree): imposed p(t) shapes, coordinate + noise model + deconvolution + window unmixing
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
SET=${SET:-r22}; O=$K_/landscape/e1b_s30_$SET; mkdir -p $O; USE=r22_; FAM=0; [ $SET = fam ] && { USE=; FAM=1; }
for sh in ${SHAPES:-sim dw2 dw4 broad shoulder flat}; do
  $PY -u $T/land_e1b.py --scores $K_/lsamp/s30a_ps2 --use "$USE" --family $FAM --sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 \
     --poses $K_/lsamp/s30a/poses_combined.txt --gt $T/simswing30/gt_particles.txt --codes swing:200:230 --basis $K_/landscape/e0_s30/truth/leaves.txt \
     --shape $sh --bstep ${BSTEP:-3} --nw ${NW:-10} --out $O/$sh 2>&1 | grep -v "Warning\|warn("
done
echo "[$(date '+%F %T')] land_e1b $SET DONE"
