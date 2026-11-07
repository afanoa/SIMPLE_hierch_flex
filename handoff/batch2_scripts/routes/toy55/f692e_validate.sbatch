#!/bin/bash
#SBATCH -J f692e_validate
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# validate the f692e states on the independent real halves: f692e leaves + closed sub-states as candidates (identity pairing),
# the cross-half pass of splitbench6 (merge test, per-shell blend pruning, cross-half EM, re-test), |CTF|.
# env PROTECT (64:65 = the closed sub-states, or empty), XS (output suffix), MINC (1500)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; E=$K/f692e; MODE=rnasep692
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; T=$K/dtree_${MODE}_f692e; O=$T/${XS:?}; mkdir -p $O; echo "f692e validation -> $O (protect '${PROTECT:-}')"
[ -s $T/cand_B.txt ] || $PY $W/f692e_candidates.py $E/consolidated_12_fullw_region/states.txt $E/consolidated_12_fullw_zoom/states.txt $W/sim$MODE $T/cand ${MINC:-1500} 2>&1 | grep -v Warning
$PY -u $W/splitbench6.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $T/cand_A.txt,$T/cand_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $O --cache $W/benchcache_$MODE --ctfabs 1 --pair identity --protect "${PROTECT:-}" --xem ${XEM:-6} --em_lim 20 --em_min 1000 ${EXTRA:-} 2>&1 | grep --line-buffered -v "Warning\|^ *warn"
for H in A B; do $PY $W/proxy_in_labels.py $K/f692c/fam1/zoomregionmotion_1_s05/level1/states_clusters.txt 1:3 $W/sim${MODE}_$H/meta_rows.txt $O/children_$H.txt "final states $H vs the f692c closed sub-states" 2>&1 | grep -v Warning; done
