#!/bin/bash
#SBATCH -J klem_xtree
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# cross-tree merge test (env MODE): the leaf states of the two independent half trees synth_<MODE>_A/_B are matched and
# merged by the cross-tree difference FSC (merge_cross_tree.py); then every half's merged states are consolidated
# (all argmax-cell members, no cap) and scored against the truth of that half
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}
RA=$A/rnasep_klem/synth_${MODE}_A; RB=$A/rnasep_klem/synth_${MODE}_B; export LEAFGLOB="region*" MASK=region UNCAP=1 MINP=${MINP:-500} NPROC=16 SPECMATCH=${SPECMATCH:-0} DOCK=${DOCK:-0} NUIS=${NUIS:-0} NULLCAL=${NULLCAL:-} FAMS=${FAMS:-?} ARCH=${ARCH:-0} THRMIN=${THRMIN:-0}; XS=${XS:-}   # XS: output suffix of a variant
echo "== (a) cross-tree merge of $MODE"; $PY $W/merge_cross_tree.py $RA $RB $W/maps160/toy_consensus.mrc $RA/cross_merge$XS 2>&1 | grep -v Warning
for H in A B; do R=$A/rnasep_klem/synth_${MODE}_$H; C=$R/consolidated_xtree$XS; NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/fam1/meta_fb.npz')['flags']))")
  echo "== (b) half $H: consolidation of the cross-tree states ($NROWS rows)"; $PY $W/consolidate_states.py $R $RA/cross_merge${XS}_$H.npz 12 $NROWS $C 2>&1 | grep -v Warning | grep -a "^cut\|dropped\|Traceback\|Error"
  echo "== (c) half $H: truth"; $PY $W/consolidated_purity.py $C $W/sim${MODE}_$H/gt_particles.txt 2>&1 | grep -v Warning; done
echo "XTREE DONE $MODE$XS"
