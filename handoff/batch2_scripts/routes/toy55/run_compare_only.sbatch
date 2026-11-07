#!/bin/bash
#SBATCH -J toy10_klem_compare
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; D=$W/lat/klem_dump20k
$PY $W/compare_estep.py $D $D/fortran_test/estep_ph0.bin 2>&1 | grep -v Warning
