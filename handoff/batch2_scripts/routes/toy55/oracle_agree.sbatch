#!/bin/bash
#SBATCH -J synth_oagree
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# oracle_agree.py over every oracle board of a set (env MODE, default arms3)
W=/mnt/beegfs/elmlund/afan/routes/toy55
/mnt/nasapps/production/python/3.12/bin/python3 $W/oracle_agree.py $W/oracle_${MODE:-arms3}/*/ 2>&1 | grep -v Warning
