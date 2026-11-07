#!/bin/bash
#SBATCH -J klem_pooled_dd
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# data-derived merge (merge_data_derived.py) of the leaves selected by LEAFGLOB (region* = level-1 leaves, zoomregion* = zoomed
# sub-states), floor FLOOR A, then consolidation (every merged state reconstructed from all its particles) and readouts
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan/rnasep_klem; export RTAG=${RTAG:-f692c} LEAFGLOB=${LEAFGLOB:-region*} FLOOR=${FLOOR:-12}
TAG=${TAG:-dd}; CUTV=$FLOOR; [ "$FLOOR" = auto ] && CUTV=12; export TREE=$A/$RTAG/merge_${TAG}_${FLOOR}.npz CTAG=_$TAG CUT=$CUTV MINP=${MINP:-1000}
echo "== data-derived merge of $LEAFGLOB leaves of $RTAG, floor $FLOOR A"; $PY $W/merge_data_derived.py $A/$RTAG $W/maps160/toy_consensus.mrc ${TREE%.npz} 2>&1 | grep -v Warning; [ -s $TREE ] || { echo 'merge FAILED'; exit 1; }
bash $W/consolidate.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^\|ValueError\|allp = "
unset LEAFGLOB; export CDIR=consolidated_${CUT}${CTAG}; bash $W/pooled_analysis.sbatch
