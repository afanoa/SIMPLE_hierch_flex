#!/bin/bash
#SBATCH -J pose_truth
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# synthetic pose accuracy vs truth for one tree (env MODE TAG FIN): consensus vs the merged leaf poses, per half
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
O=$K/dtree_${MODE}_$TAG/${FIN:-final_xem}
for h in A B; do echo "== $TAG half $h"
  $PY $W/pose_truth_check.py $W/sim${MODE}_$h/all.txt $W/sim${MODE}_$h/gt_particles.txt $O/leaves_$h.txt consensus=$W/lat_${MODE}_$h/ptcl3D_r0.txt merged_refined=$O/poses_$h.txt $EXTRA_$h 2>&1 | grep -v Warning
done
