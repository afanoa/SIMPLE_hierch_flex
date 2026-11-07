#!/bin/bash
#SBATCH -J klem_track
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# track_node.py on toy trees against simsites truth: where REF's node NODE went in each tree; env REF NODE TREES 'name:tree_dir ..', OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
ST=$K/factorial/junk_test/simsites_truth_jt.txt; OUT=${OUT:?}; mkdir -p $(dirname $OUT)
for e in ${TREES:?}; do IFS=: read nm dir <<< "$e"; echo "######## $nm"; $PY $T/track_node.py ${REF:?} ${NODE:?} $dir $ST 2>&1 | grep -v Warning; done | tee $OUT
