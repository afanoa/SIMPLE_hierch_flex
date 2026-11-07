#!/bin/bash
#SBATCH -J toy10_modelclass2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Model-class test rerun with weaker field regularisation (the first pass with sigma 8/14 px barely moved: 1-2 A)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W
MAPS=$(ls $W/maps160/docked/conf*_m160.mrc | sort | tr "\n" , | sed 's/,$//'); NAMES=$(ls $W/maps160/docked/conf*_m160.mrc | sort | xargs -n1 basename | cut -d_ -f1 | tr "\n" , | sed 's/,$//')
MAPS2="$MAPS,$(ls $W/maps160/variants/conf*_m160.mrc | sort | tr "\n" , | sed 's/,$//')"; NAMES2="$NAMES,$(ls $W/maps160/variants/conf*_m160.mrc | sort | xargs -n1 basename | cut -d_ -f1 | tr "\n" , | sed 's/,$//')"
for SG in 3 6; do mkdir -p $W/model_class_s$SG
  $PY $W/model_class.py $W/model_class_s$SG $W/maps160/toy_consensus.mrc $MAPS --names $NAMES --nproc 8 --sigma $SG 2>&1 | tee $W/model_class_s$SG/log.txt | grep -v "^   linear LOO" | cut -c1-400; done
mkdir -p $W/model_class_c3
$PY $W/model_class.py $W/model_class_c3 $W/maps160/toy_consensus.mrc $MAPS2 --names $NAMES2 --nproc 8 --sigma 3 2>&1 | tee $W/model_class_c3/log.txt | grep -v "^   linear LOO" | cut -c1-400
echo DONE
