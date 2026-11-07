#!/bin/bash
#SBATCH -J klem_wholefac
#SBATCH -p norm
#SBATCH --qos=norm
#SBATCH -c 64
#SBATCH --mem=150G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# a global motion (e.g. breathing of the whole periphery) is invisible to local region boards: add one whole-particle board per family
# of a fast run (no focus mask, region_c0) as a candidate factor, then factor_sample.py with the axis readout. env TAG, RT (ff2_<TAG>), CODES
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; W=$A/routes/toy55_f692v2; K=$A/rnasep_klem
TAG=${TAG:?}; R=$K/${RT:-ff2_$TAG}; NF=$(cat $R/expand/nfam.txt); export FB=$K/f692e/klem2_bin CTFABS=0 PSCALE=0 SERVER_LIST= NPARTS=
for F in $(seq 1 $NF); do D=$R/fam$F; S=$R/expand/family${F}_subset.txt; mkdir -p $D/region_c0
  [ -s $D/region_c0/board/klem_latent.txt ] || ( NTHR=21 N1=7 N2=7 NFREE=0 META=$R/meta_toy.npz bash $W/fboard.sh $D/region_c0/board $R/whole/p.simple $S "12 10 10 8 8" "1.2 1.0 1.0 0.8 0.8" > $D/region_c0/board.out 2>&1 ) &
done; wait; echo "whole-particle family boards: $(ls $R/fam*/region_c0/board/klem_latent.txt 2>/dev/null | wc -l) of $NF"
$PY $T/factor_sample.py $R $R/factors_whole --fams $(seq -s, 1 $NF) --glob 'region_c*' --nbin 8 --budget 1000 --minend 300 --gt $T/sim$TAG/gt_particles.txt --codes ${CODES:?} 2>&1 | grep -v Warning | grep 'family\|factors at\|factor [0-9]\|rank correlation\|region_c0 '
