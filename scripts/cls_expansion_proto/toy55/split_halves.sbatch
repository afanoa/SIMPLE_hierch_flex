#!/bin/bash
#SBATCH -J split_halves
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/split_halves.%j.out
# independent halves of a synthetic set (env MODE, NP): sim<MODE>_A/_B, then poses (run_toy55_scale, no simulation)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; NP=${NP:-10}; cd $W
$PY $W/split_halves.py $W/sim$MODE $W/sim${MODE}_A $W/sim${MODE}_B || exit 1
for H in A B; do sbatch -c ${MC:-48} --export=ALL,MODE=${MODE}_$H,N=1,NP=$NP run_toy55_scale.sbatch; done
echo "SPLIT DONE $MODE"
