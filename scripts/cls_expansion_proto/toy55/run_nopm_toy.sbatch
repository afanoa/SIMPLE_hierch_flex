#!/bin/bash
#SBATCH -J toy10_nopm
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# A/B of pose marginalisation inside the final recipe: the balanced realigned tree's children and level 3 WITHOUT
# --pose_marg (one E-step instead of seven), against 0.64 (children) and 0.71 (level 3) with it.
set -uo pipefail
W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') nopm] $*"; }
log "children without pose marginalisation"; env RTAG=_realign1_bal05_nopm LOOP=0 BAL=0.5 PM=0 bash $W/run_tree_realign_toy.sbatch 2>&1 | grep "owning at least\|MERGED\|best-state\|Error\|Traceback\|DONE" | sed 's/^/[children] /'
log "level 3 without pose marginalisation"; env CB=$W/lat/klem_expand_realign1_bal05_nopm_fam TABLE=$W/lat/realign_realign1/ptcl3D_realigned.txt BAL=0.5 PMARG=0 L3=_l3b_nopm bash $W/run_tree_level3b_toy.sbatch 2>&1 | grep "states: cross-half\|MERGED\|best-state\|Error\|Traceback\|DONE" | sed 's/^/[level3] /'
log DONE
