#!/bin/bash
#SBATCH -J klem_maps_arm
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# final-state maps from the bench geometry (both halves averaged, full box, frame checked) + closed-arm readout + state graph and
# membership UMAP (env MODE, TAG, FIN, POSEA/POSEB (default the lat tables), CONS, CLOSED, EXTRA e.g. "--ctfabs 1", NOGT)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}; FIN=${FIN:-final_xem}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; O=$K/dtree_${MODE}_$TAG/$FIN; P=$K/graph_${MODE}_${TAG}_$FIN
TR=""; [ "${MODE#55}" != "$MODE" ] && TR="--truth $W/truthtree_55_200k.npz"; GT=""; [ -s $W/sim${MODE}_A/gt_particles.txt ] && [ "${NOGT:-0}" = 0 ] && GT="--gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt"
$PY $W/state_graph.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses ${POSEA:-$W/lat_${MODE}_A/ptcl3D_r0.txt},${POSEB:-$W/lat_${MODE}_B/ptcl3D_r0.txt} --labels $O/children_A.txt,$O/children_B.txt \
  --cache $W/benchcache_$MODE --png $P.png --write_maps $O/maps_bench --cons ${CONS:?} $GT $TR ${EXTRA:-} 2>&1 | grep --line-buffered -v Warning
echo "== closed-arm readout of the final states"; $PY $W/arm_readout.py $CONS ${CLOSED:?} $O/arm_readout_bench.txt $O/maps_bench 2>&1 | grep -v Warning | sed "s#$O/##"
$PY $W/umap_states.py $P.npz ${P}_umap.png 2>&1 | grep -v "Warning\|warn(" | tail -1
