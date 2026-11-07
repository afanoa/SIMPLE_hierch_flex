#!/bin/bash
#SBATCH -J land_confs
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_confs.%j.out
# deposited conformers, holo conformers, rigid-swing maps and f692st leaves: arm readouts and positions on the rigid swing
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=$K_/landscape/confs; mkdir -p $O
SW=$T/maps160/swing60
[ -s $SW/truth_swing.txt ] || $PY $T/make_swing_maps.py $T/maps160/ah_norm/conf18_emd74209_m160.mrc $T/maps160/variants/domain_mask.mrc $SW \
  --open $(ls $T/maps160/ah_norm/conf27_*_m160.mrc $T/maps160/ah_norm/conf28_*_m160.mrc $T/maps160/ah_norm/conf48_*_m160.mrc $T/maps160/ah_norm/conf49_*_m160.mrc | tr '\n' ',' | sed 's/,$//') \
  --travel 60 --nsteps 31 --id0 400 2>&1 | grep --line-buffered -v Warning
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_all.txt $T/maps160/docked_all $T/maps160/holo_docked $SW $K_/rflex_f692st/leafmaps 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | tail -3
for sw in swing30 swing60; do $PY $T/land_swingfit.py $T/maps160/$sw $O/swingfit_$sw.txt $T/maps160/docked_all $T/maps160/holo_docked $K_/rflex_f692st/leafmaps $T/maps160/swing60 2>&1 | grep --line-buffered -v Warning | head -2; done
echo "[$(date '+%F %T')] land_confs DONE"
