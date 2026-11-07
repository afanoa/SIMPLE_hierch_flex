#!/bin/bash
#SBATCH -J avgs_partial
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# per-Fab readout of finished sub-class E-step runs (env MODE, RUNS = run dirs under avgs_<MODE> or bags_<MODE>)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:-armsigm2}
for R in ${RUNS:?}; do echo "== $(basename $R)"; $PY $W/factorised_readout.py $R $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning | grep -a "per arm"; done; echo PARTIAL DONE
