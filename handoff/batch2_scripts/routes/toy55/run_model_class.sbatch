#!/bin/bash
#SBATCH -J toy10_modelclass
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Stream A of the latent-variable night: model-class test on the 10 docked truth conformers (linear span vs smooth
# deformation vs deformation + occupancy, per band), then the two compositional variants (mobile domain deleted,
# extra density added), the same test on the 12 maps, and the launch of the compositional toy simulation (toy10c).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W
log(){ echo "[$(date '+%F %T') modelclass] $*"; }
n=0; while [ $n -lt 240 ] && [ ! -s $W/maps160/toy_consensus.mrc ]; do sleep 60; n=$((n+1)); done
[ -s $W/maps160/toy_consensus.mrc ] || { log "no toy consensus"; exit 1; }
MAPS=$(ls $W/maps160/docked/conf*_m160.mrc | sort | tr "\n" , | sed 's/,$//'); NAMES=$(ls $W/maps160/docked/conf*_m160.mrc | sort | xargs -n1 basename | cut -d_ -f1 | tr "\n" , | sed 's/,$//')
log "START: 10 conformers ($NAMES)"
mkdir -p $W/model_class
$PY $W/model_class.py $W/model_class $W/maps160/toy_consensus.mrc $MAPS --names $NAMES --nproc 8 --sigma 8 2>&1 | tee $W/model_class/log.txt | grep -v "^   linear LOO" | cut -c1-400
# stiffer field (sigma 14 px = 20 A): closer to a rigid-body/domain model
mkdir -p $W/model_class_s14
$PY $W/model_class.py $W/model_class_s14 $W/maps160/toy_consensus.mrc $MAPS --names $NAMES --nproc 8 --sigma 14 2>&1 | tee $W/model_class_s14/log.txt | grep -A9 "^SUMMARY\|^flow-field" | cut -c1-300
log "compositional variants"
$PY $W/make_variants.py $W/rnasep_conformers10.csv $W/maps160/docked $W/maps160/toy_consensus.mrc $W/maps160/variants $W/rnasep_conformers10c.csv 2>&1 | tee $W/maps160/variants/log.txt
MAPS2="$MAPS,$(ls $W/maps160/variants/conf*_m160.mrc | sort | tr "\n" , | sed 's/,$//')"; NAMES2="$NAMES,$(ls $W/maps160/variants/conf*_m160.mrc | sort | xargs -n1 basename | cut -d_ -f1 | tr "\n" , | sed 's/,$//')"
log "launching the compositional toy simulation (toy10c)"
sbatch $W/run_toy10c_sim.sbatch
mkdir -p $W/model_class_c
$PY $W/model_class.py $W/model_class_c $W/maps160/toy_consensus.mrc $MAPS2 --names $NAMES2 --nproc 8 --sigma 8 2>&1 | tee $W/model_class_c/log.txt | grep -v "^   linear LOO" | cut -c1-400
log "DONE stream A"
