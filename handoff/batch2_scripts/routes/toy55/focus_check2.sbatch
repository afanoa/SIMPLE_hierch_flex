#!/bin/bash
#SBATCH -J focus_check
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the pcg_focus invariant (focus_check2.py) on a finished fam_classify_focus run: env RD (its refine3D dir), MASK, K (8)
set -euo pipefail
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/focus_check2.py ${RD:?} ${MASK:?} ${K:-8}
echo DONE
