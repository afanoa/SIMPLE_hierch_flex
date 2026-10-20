#!/bin/bash
#SBATCH -J klem_belong
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# do the leaves of a data-derived tree belong to their parents? cross-half map and particle checks (env MODE, TAG)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; O=$K/dtree_${MODE}_$TAG/final; echo "tree $MODE $TAG"
TR=""; [ "${MODE#55}" != "$MODE" ] && TR="--truth $W/truthtree_55_200k.npz"
$PY $W/tree_belong.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --leaves $O/leaves_A.txt,$O/leaves_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --cache $W/benchcache_$MODE $TR ${EXTRA:-} 2>&1 | grep -v Warning
