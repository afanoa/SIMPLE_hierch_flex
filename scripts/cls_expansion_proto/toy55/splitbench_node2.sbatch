#!/bin/bash
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# one node's split test for dtree.sbatch (BACKEND=bench): args CAND_A CAND_B OUTDIR, env MODE, NSHUF (20), CROP (48)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; mkdir -p $3
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
$PY $W/splitbench2.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $W/lat_${MODE}_A/ptcl3D_r0.txt,$W/lat_${MODE}_B/ptcl3D_r0.txt --cand $1,$2 \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $3 --nshuf ${NSHUF:-20} --crop ${CROP:-48} --cache $W/benchcache_$MODE 2>&1 | grep -v Warning > $3/bench.log
