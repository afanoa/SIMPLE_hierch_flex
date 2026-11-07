#!/bin/bash
#SBATCH -J klem_bench
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# fast in-memory split test of one data-tree node (env MODE, TAG (dtree2), ND (r), NSHUF (20), CROP (48), EXTRA)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:-dtree2}; ND=${ND:-r}
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
C=fam0/regioncand/level1/states_full.txt; CA=$K/synth_${MODE}_A/$TAG/$ND/$C; CB=$K/synth_${MODE}_B/$TAG/$ND/$C
echo "== splitbench $MODE $TAG node $ND"; t0=$(date +%s)
$PY $W/splitbench2_shell.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $CA,$CB \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $K/splitbench_${MODE}_${TAG}_$ND --nshuf ${NSHUF:-20} --crop ${CROP:-48} --cache $W/benchcache_$MODE ${EXTRA:-} 2>&1 | grep -v Warning
echo "[$(date '+%F %T') bench] done in $(( $(date +%s) - t0 )) s"
