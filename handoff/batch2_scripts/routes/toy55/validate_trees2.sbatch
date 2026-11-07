#!/bin/bash
#SBATCH -J klem_validate2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# validation of port toy trees against truth (validate_tree2.py), one SUMMARY line per tree; env TREES 'name:kind:tree_dir[:table] ..'
# (kind sites = simsites truth, arms = IgM truth), OUT (summary file)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
ST=$K/factorial/junk_test/simsites_truth_jt.txt; AT=$K/factorial/armsigm1_truth/arms_truth.txt; OUT=${OUT:-$K/factorial/validation/summary.txt}; mkdir -p $(dirname $OUT)
for e in ${TREES:?}; do IFS=: read nm kind dir tab <<< "$e"; tr=$ST; [ "$kind" = arms ] && tr=$AT
  $PY $T/validate_tree2.py $nm $dir $tr $kind ${tab:+--table $tab} 2>&1 | grep -v Warning; done | tee $OUT.full | grep "^SUMMARY" > $OUT
echo "######## summary"; cat $OUT
