#!/bin/bash
#SBATCH -J klem_splitsup
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# seed recurrence of splits (split_support.py) for pairs of trees on the same particles; env PAIRS 'kind:name1:tree1:name2:tree2 ..'
# (kind sites = simsites truth, conf = CAHRA 13991 STAR truth), OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
ST=$K/factorial/junk_test/simsites_truth_jt.txt; GT=/mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star; OUT=${OUT:?}
for e in ${PAIRS:?}; do IFS=: read kind n1 d1 n2 d2 <<< "$e"; tr=$ST; [ "$kind" = conf ] && tr=$GT
  echo "######## $n1 vs $n2"; $PY $T/split_support.py $kind $tr $n1 $d1 $n2 $d2 2>&1 | grep -v Warning; done | tee $OUT
