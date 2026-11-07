#!/bin/bash
#SBATCH -J klem_recount
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape recount of a finished run (env RTAG) after late regions: the chain's step (5) over every leaf board
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}
DIRS=""; for d in $R/fam?/klem_fb $R/fam?/region*/board; do [ -s $d/klem_halves.npz ] && DIRS="$DIRS $d"; done
for d in $DIRS; do [ -s $d/cell_matches.txt ] || $PY $W/deposited_match.py $d/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $d/cell_matches.txt > /dev/null 2>&1; done
echo "== LANDSCAPE over the leaves ($RTAG, $(ls $R/fam?/region*/DONE 2>/dev/null | wc -l) regions DONE): $DIRS"; $PY $W/count_states.py $DIRS 2>&1 | grep -v Warning | tail -20
