#!/bin/bash
#SBATCH -J klem_count_states_real
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Deposited matches saved for every real leaf board, then the distinct-state count with deposited identification per
# cluster: family 1 (feedback board + level-3 grandchildren) and, when present, family 2 (feedback board + grandchildren).
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan; E=$A/rnasep_klem/expand1
DIRS="$E/fam1_pipeline/klem_fb $E/fam1_l3/sub1 $E/fam1_l3/sub2"; for d in $E/fam2_pipeline/klem_fb $E/fam2_l3/sub1 $E/fam2_l3/sub2; do [ -s $d/klem_halves.npz ] && DIRS="$DIRS $d"; done
for d in $DIRS; do [ -s $d/cell_matches.txt ] || $PY $W/deposited_match.py $d/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $d/cell_matches.txt > /dev/null 2>&1; done
echo "== REAL leaves with deposited identification: $DIRS"; $PY $W/count_states.py $DIRS 2>&1 | grep -v Warning
