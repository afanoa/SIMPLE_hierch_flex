#!/bin/bash
#SBATCH -J klem_zoom_regmerge
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# region-local merge of the zoomed sub-states (deviations masked by the region of each leaf), then consolidation at CUT
# (union of particles per group reconstructed), arm readout and strict recovery of 55
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692c}; export LEAFGLOB='zoomregion*' REGMASK=1
echo "== region-local map-similarity tree of the zoomed sub-states"; $PY $W/merge_states.py $R $W/maps160/toy_consensus.mrc $R/merge_tree_zoom_reg 2>&1 | grep -v Warning | grep -E "similarity|states|region mask" | tail -20
echo REGMERGE DONE
