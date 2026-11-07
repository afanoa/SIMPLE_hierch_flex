#!/bin/bash
#SBATCH -J klem_rflex12test
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# node test of rflex12 (borrowed stage B): a NEW run RT whose root node starts from copies of a finished run's root node up to stage A (FROM:
# n.simple, mean, flex rows and states, rsmask, focused flex, stage A), so only stage B (BORROW), the split test and the decision run (MAXD 0:
# the supported groups are leaves); then the closed-arm readout of the stage-B states (node_readout2.sbatch) and of the leaves (rflex12's own).
# env FROM (finished run tag), RT (new tag), BORROW (auto | none), SRC, ROOTROWS, MSKDIAM, KFLOOR; rflex12's other keys as usual
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; S=$K_/rflex_${FROM:?}/node_r; R=$K_/rflex_${RT:?}; D=$R/node_r
[ -e $R ] && [ ! -e $R/seeded ] && { echo "$R exists and was not seeded by this script: refusing"; exit 1; }
if [ ! -e $R/seeded ]; then mkdir -p $D && cp -a $S/n.simple $S/flex_rows.txt $S/rsmask.mrc $S/rsmask_variance.mrc $D/ && cp -a $S/mean $S/1_flex_pca $S/focus $S/fstageA $D/ && touch $R/seeded || { echo "seeding FAILED"; exit 1; }; fi
echo "[$(date '+%F %T')] $RT seeded from $FROM (stage A), BORROW ${BORROW:-auto}"
MAXD=${MAXD:-0} bash $T/rflex12.sbatch
RD=$(ls -td $D/fstageB/*_refine3D | head -1) LAB=$D/refined_f.txt OUT=$R/stageB_readout VIEW=$K_/view_rflex/${RT}_r_fstageB NAME="$RT stage B (BORROW ${BORROW:-auto})" bash $T/node_readout2.sbatch
echo "[$(date '+%F %T')] $RT DONE"
