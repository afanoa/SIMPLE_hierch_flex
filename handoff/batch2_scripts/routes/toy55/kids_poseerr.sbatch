#!/bin/bash
#SBATCH -J kids_poseerr
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Trap check from the split-decision session: do the cheap PC1 children sort by the particles' pose ERROR (truth simswing30/all.txt vs
# the parent's poses)? pose_err_groups.py per child: mean 3-D shift-error and rotation-error vectors vs their standard errors.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=4
for u in r22 r11 r12 r21; do echo "===== s30a $u cheap children"; D=$K_/cheapkids_s30a_$u
  grep "e1=" $D/ptcl3D_parent.txt > $D/cons_oris.txt; $PY -u $T/pose_err_groups.py $T/simswing30/all.txt $D/cons_oris.txt $D/children.txt 2>&1 | grep --line-buffered -v Warning; done
echo "[$(date '+%F %T') kids_poseerr] DONE"
