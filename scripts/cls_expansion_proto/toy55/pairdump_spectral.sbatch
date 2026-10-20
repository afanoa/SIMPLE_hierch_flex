#!/bin/bash
#SBATCH -J klem_pspec
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/pairdump_spectral.py /mnt/beegfs/elmlund/afan/rnasep_klem/pairdump_null_50k_r_regionshuf1_max.npz 2>&1 | grep -v Warning
