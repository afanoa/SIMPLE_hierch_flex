#!/bin/bash
#SBATCH -J focus_mask
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the swept-domain focus mask of ff3_swing30hs (swing_mask.py) into rnasep_klem/ff3_swing30hs/focus/swing_swept.mrc
set -euo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; O=$A/rnasep_klem/ff3_swing30hs/focus
mkdir -p $O
/mnt/nasapps/production/python/3.12/bin/python3 $T/swing_mask.py $T/maps160/swing30hs $O/swing_swept.mrc \
    $A/rnasep_klem/ff2_swing30hs/cand/region_5.mrc
echo DONE
