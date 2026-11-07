#!/bin/bash
#SBATCH -J klem_simmat
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# leaf x leaf cross-half similarity of a data-derived tree, ordered by the tree (env MODE, TAG)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; O=$K/dtree_${MODE}_$TAG/final; echo "tree $MODE $TAG"
TR=""; [ "${MODE#55}" != "$MODE" ] && TR="--truth $W/truthtree_55_200k.npz"
$PY $W/tree_simmat.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --leaves $O/leaves_A.txt,$O/leaves_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --cache $W/benchcache_$MODE --png $K/simmat_${MODE}_$TAG${XS:-}.png ${EXTRA:-} $TR 2>&1 | grep --line-buffered -v Warning
