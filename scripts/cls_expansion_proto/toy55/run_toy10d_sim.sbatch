#!/bin/bash
#SBATCH -J toy10d_sim
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# toy10d: the 10 deposited conformers + conf63 (most populous conformer with 35 % of its mass cut away: the assembly-intermediate
# regime) at 15 % of the conformer particles; then the standard toy simulation chain (run_toy10c_sim.sbatch with simd/projd).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
mkdir -p $W/maps160/variants_big
$PY $W/make_big_variant.py $W/rnasep_conformers10.csv $W/maps160/docked $W/maps160/variants_big $W/rnasep_conformers10d.csv --frac ${FRAC:-0.35} --images ${IMAGES:-186000} || exit 1
mkdir -p $W/simd; for c in sim simc; do [ -d $W/simd/calib ] || { [ -d $W/$c/calib ] && cp -r $W/$c/calib $W/simd/calib; }; done   # the SNR calibration uses the toy consensus only
export SIMDIR=simd PROJ=projd CSV=rnasep_conformers10d.csv TAG=d EXTRA_MAPS="$W/maps160/variants_big/conf63_*_m160.mrc" N=${N:-200000}
bash $W/run_toy10c_sim.sbatch
