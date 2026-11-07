#!/bin/bash
#SBATCH -J klem_tuples
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/tuple_landscape.py /mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:-f692e} 932496 --minp ${MINP:-1000} 2>&1 | grep -v Warning; echo "TUPLES DONE"
