#!/bin/bash
#SBATCH -J klem_umap
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# state graph (saves the cross-half memberships) then the per-particle UMAP of those memberships (env MODE, TAG, FIN, EXTRA, NOGT)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}; FIN=${FIN:-final_xem}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; O=$K/dtree_${MODE}_$TAG/$FIN; P=$K/graph_${MODE}_${TAG}_$FIN
TR=""; [ "${MODE#55}" != "$MODE" ] && TR="--truth $W/truthtree_55_200k.npz"; GT=""; [ -s $W/sim${MODE}_A/gt_particles.txt ] && [ "${NOGT:-0}" = 0 ] && GT="--gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt"
$PY $W/state_graph.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --labels $O/children_A.txt,$O/children_B.txt \
  --cache $W/benchcache_$MODE --png $P.png $GT $TR ${EXTRA:-} 2>&1 | grep --line-buffered -v Warning
$PY $W/umap_states.py $P.npz ${P}_umap.png 2>&1 | grep -v Warning | tail -3
