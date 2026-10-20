#!/bin/bash
#SBATCH -J klem_bench
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-shell archetype bench on a node of a CONFS-cut tree (env MODE, TAG, ND)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
C=fam0/regioncand/level1/states_full.txt; echo "== splitbench (archmode shell) $MODE $TAG node $ND"
$PY $W/splitbench2_shell.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $K/synth_${MODE}_A/$TAG/$ND/$C,$K/synth_${MODE}_B/$TAG/$ND/$C --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $K/splitbench_${MODE}_${TAG}_${ND}_shell --cache $W/benchcache_$MODE --archmode shell 2>&1 | grep -v Warning
for H in A B; do echo " half $H clades:"; $PY $W/clade_compose.py $K/splitbench_${MODE}_${TAG}_${ND}_shell/children_$H.txt $W/sim${MODE}_$H/gt_particles.txt $W/truthtree_55_200k.npz 0.45 2>&1 | grep -v "Warning\|reference clades"; done
