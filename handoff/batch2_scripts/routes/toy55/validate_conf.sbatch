#!/bin/bash
#SBATCH -J klem_validate_conf
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# validation of trees on a two-conformation STAR truth (validate_conf.py), one SUMMARY line per tree; env TREES 'name:tree_dir[:replicate_dir] ..',
# GT (truth STAR, default CAHRA 13991 source.star), OUT (summary file)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
GT=${GT:-/mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star}; OUT=${OUT:-$K/factorial/validation/cahra.txt}; mkdir -p $(dirname $OUT)
for e in ${TREES:?}; do IFS=: read nm dir rep <<< "$e"
  $PY $T/validate_conf.py $nm $dir $GT ${rep:+--rep $rep} 2>&1 | grep -v Warning; done | tee $OUT.full | grep "^SUMMARY" > $OUT
echo "######## summary"; cat $OUT
