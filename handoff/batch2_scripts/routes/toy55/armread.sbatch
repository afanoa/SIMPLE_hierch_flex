#!/bin/bash
#SBATCH -J armread
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# per-arm truth readout of synthetic klem trees (env MODE and RUNS = run tags under rnasep_klem, space separated)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}
for t in ${RUNS:?}; do echo "=== $t"; $PY $W/arm_angle_readout.py $A/rnasep_klem/$t $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning; done
echo ARMREAD DONE
