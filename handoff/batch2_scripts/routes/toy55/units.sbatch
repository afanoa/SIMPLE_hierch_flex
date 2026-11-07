#!/bin/bash
#SBATCH -J klem_units
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merge-up units scored like leaves: units_prep.py builds each unit layout, validate_tree3.py scores it against simsites truth;
# env TREES 'name:tree_dir[:state_labels] ..' (state_labels = E-step states mapped through the tree's leaf -> unit map), OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
ST=$K/factorial/junk_test/simsites_truth_jt.txt; OUT=${OUT:?}; W=$(dirname $OUT)/units; mkdir -p $W
for e in ${TREES:?}; do IFS=: read nm dir lab <<< "$e"
  $PY $T/units_prep.py $dir $W/$nm ${lab:+--labels $lab} 2>&1 | grep -v Warning
  $PY $T/validate_tree3.py $nm $W/$nm $ST sites 2>&1 | grep -v Warning | grep -v "   split r "; done | tee $OUT.full | grep "^SUMMARY" > $OUT
echo "######## summary"; cat $OUT
