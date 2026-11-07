#!/bin/bash
#SBATCH -J klem_rflex13test2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rflex13 node test with its own FALLOFF mask, borrowed vs plain stage B on identical inputs. Phase 1: new run RTN (BORROW none) seeded with
# FROM's root node only up to the plain flex and the cheap children (n.simple, mean, flex rows and states, cheap/), so rflex13 builds the falloff
# mask (arm_mask2.py), the focused flex, stage A, stage B, the focus-mask split test and the leaves (MAXD 0). Phase 2: new run RTB (BORROW auto)
# seeded from RTN's root up to stage A (same mask, focused flex, stage A): stage B borrowed + split test. Readouts: stage-B states
# (node_readout2.sbatch: closed-arm, old f692r closed 15+16) and rflex13's leaf readout for both. env FROM, RTN, RTB, SRC, ROOTROWS, MSKDIAM, KFLOOR
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; S=$K_/rflex_${FROM:?}/node_r; RN=$K_/rflex_${RTN:?}; RB=$K_/rflex_${RTB:?}
for R in $RN $RB; do [ -e $R ] && [ ! -e $R/seeded ] && { echo "$R exists and was not seeded by this script: refusing"; exit 1; }; done
if [ ! -e $RN/seeded ]; then mkdir -p $RN/node_r && cp -a $S/n.simple $S/flex_rows.txt $RN/node_r/ && cp -a $S/mean $S/1_flex_pca $S/cheap $RN/node_r/ && touch $RN/seeded || { echo "phase 1 seeding FAILED"; exit 1; }; fi
echo "[$(date '+%F %T')] phase 1: $RTN seeded from $FROM (plain flex + cheap children), BORROW none"
RT=$RTN BORROW=none MAXD=${MAXD:-0} bash $T/rflex13.sbatch
[ -e $RN/node_r/node.done ] || { echo "phase 1 did not finish"; exit 1; }
RD=$(ls -td $RN/node_r/fstageB/*_refine3D | head -1) LAB=$RN/node_r/refined_f.txt OUT=$RN/stageB_readout VIEW=$K_/view_rflex/${RTN}_r_fstageB NAME="$RTN stage B (BORROW none, falloff mask)" bash $T/node_readout2.sbatch
if [ ! -e $RB/seeded ]; then D=$RN/node_r; mkdir -p $RB/node_r && cp -a $D/n.simple $D/flex_rows.txt $D/rsmask.mrc $RB/node_r/ && { [ -e $D/rsmask_variance.mrc ] && cp -a $D/rsmask_variance.mrc $RB/node_r/ || true; } \
  && cp -a $D/mean $D/1_flex_pca $D/focus $D/fstageA $RB/node_r/ && touch $RB/seeded || { echo "phase 2 seeding FAILED"; exit 1; }; fi
echo "[$(date '+%F %T')] phase 2: $RTB seeded from $RTN (stage A, falloff mask), BORROW auto"
RT=$RTB BORROW=auto MAXD=${MAXD:-0} bash $T/rflex13.sbatch
RD=$(ls -td $RB/node_r/fstageB/*_refine3D | head -1) LAB=$RB/node_r/refined_f.txt OUT=$RB/stageB_readout VIEW=$K_/view_rflex/${RTB}_r_fstageB NAME="$RTB stage B (BORROW auto, falloff mask)" bash $T/node_readout2.sbatch
echo "[$(date '+%F %T')] DONE"
