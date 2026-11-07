#!/bin/bash
#SBATCH -J klem_finalize
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# final consolidation of a data-derived tree: all leaves of both halves through the split test (pairing, merges, per-shell
# blend pruning, soft EM) -> the number of states (env MODE, TAG)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; O=$K/dtree_${MODE}_$TAG/final; mkdir -p $O
for H in A B; do $PY $W/tree_leaves.py $K/synth_${MODE}_$H $TAG $O/leaves_$H.txt; done
$PY $W/splitbench3.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $O/leaves_A.txt,$O/leaves_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $O --cache $W/benchcache_$MODE ${EXTRA:-} 2>&1 | grep -v Warning
if [ "${MODE#55}" != "$MODE" ]; then for H in A B; do echo " half $H final states vs reference clades:"; $PY $W/clade_compose.py $O/children_$H.txt $W/sim${MODE}_$H/gt_particles.txt $W/truthtree_55_200k.npz 0.45,0.40 2>&1 | grep -v Warning; done; fi
