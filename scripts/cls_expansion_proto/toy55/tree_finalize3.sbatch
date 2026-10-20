#!/bin/bash
#SBATCH -J klem_finalize3
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# finalize3 = finalize2 with splitbench5 (--ctfabs CTFABS for phase-flipped real data; NOGT=1 skips the truth readouts). final consolidation with a CROSS-HALF re-assignment: all leaves of both halves through the split test,
# cross-half EM (each half against the other half's filtered group maps), re-test -> the number of states; then the truth
# readout of the leaves, the old final pass (own-half soft EM) and this one (env MODE, TAG; XEM iterations, EMLIM band in A)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; D=$K/dtree_${MODE}_$TAG; O=$D/final_xem${XS:-}; mkdir -p $O; echo "tree $MODE $TAG: cross-half final pass, $O"
for H in A B; do $PY $W/tree_leaves.py $K/synth_${MODE}_$H $TAG $O/leaves_$H.txt; done
$PY -u $W/splitbench5.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $O/leaves_A.txt,$O/leaves_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $O --cache $W/benchcache_$MODE --xem ${XEM:-6} --em_lim ${EMLIM:-20} --em_min ${EMMIN:-1000} --ctfabs ${CTFABS:-0} ${EXTRA:-} 2>&1 | grep --line-buffered -v Warning
TR=""; [ "${MODE#55}" != "$MODE" ] && TR="$W/truthtree_55_200k.npz 0.45"; GT=$W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt
if [ "${NOGT:-0}" = 1 ]; then exit 0; fi
echo "== truth readout: tree leaves, old final pass (own-half soft EM), new final pass (cross-half EM)"
$PY $W/final_nmi.py leaves $O/leaves_A.txt,$O/leaves_B.txt $GT $TR 2>&1 | grep -v Warning
[ -s $D/final/children_A.txt ] && $PY $W/final_nmi.py old_final $D/final/children_A.txt,$D/final/children_B.txt $GT $TR 2>&1 | grep -v Warning
$PY $W/final_nmi.py new_final $O/children_A.txt,$O/children_B.txt $GT $TR 2>&1 | grep -v Warning
