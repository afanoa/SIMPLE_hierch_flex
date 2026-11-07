#!/bin/bash
#SBATCH -J klem_sbgrid
#SBATCH -p norm
#SBATCH --qos=norm
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the swing x breathing map grid once (make_swingbreath_maps.py into maps160/sbgrid), then both product simulations (indep, coupled)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; M=$W/maps160/ah_norm; cd $W
[ -s $W/maps160/sbgrid/truth_grid.txt ] || $PY $W/make_swingbreath_maps.py $M/conf18_emd74209_m160.mrc $W/maps160/variants/domain_mask.mrc $W/maps160/sbgrid \
  --open $(ls $M/conf27_*_m160.mrc $M/conf28_*_m160.mrc $M/conf48_*_m160.mrc $M/conf49_*_m160.mrc | tr '\n' ',' | sed 's/,$//') 2>&1 | grep -v Warning
[ -s $W/maps160/sbgrid/truth_grid.txt ] || { echo "grid maps FAILED"; exit 1; }
for m in indep coupled; do j=$(sbatch --parsable --export=ALL,MODE=$m run_product_sim.sbatch); echo "sim $m $j"; sbatch --parsable --export=ALL,SIMJ=$j,TAG=prod_$m,RTAG=f692p${m:0:1} prep_swing_chain.sbatch; done
