#!/bin/bash
#SBATCH --job-name=leafdiag
#SBATCH --partition=csbdevel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/beegfs/elmlund/afan/rnasep_klem/leafdiag.%j.out
module load gcc/15.1.0 2>/dev/null; W=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3
M=${M:-null_200k}; C=${C:-/mnt/beegfs/elmlund/afan/rnasep_klem/synth_$M/consolidated_12_vm}
echo "== leaf diagnostics $M ($C)"; $PY $W/leaf_diag.py $W/sim$M/all.txt $W/lat_$M/ptcl3D_r0.txt ${ST:-$C/states.txt}
