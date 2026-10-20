#!/bin/bash
#SBATCH -J rnasep_f692_docked
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Docked-node enumeration after the f692 chain (f692_docked.sh): env RTAG (f692), FAM (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692}; FAM=${FAM:-1}
export R FAM META=$R/fam$FAM/meta_fb.npz FB=$R/klem2_bin
bash $W/f692_docked.sh
