#!/bin/bash
#SBATCH -J synth_llrbias
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# column-centred readout of an existing arm_llr.npz (env NPZ)
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/llr_bias.py ${NPZ:?} 2>&1 | grep -v Warning
