#!/bin/bash
#SBATCH -J klem_rflex_post
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# post-processing of an rflex run (rflex.sbatch): the leaves as integer ids -> their bottom-up tree (state_merge.py --list) ->
# tree-weighted leaf maps against truth (tree_wiener4.py --list: own / star / tree / tree2 / emp2 ...) -> the truth readouts.
# env RT, TAG, MSKDIAM (220)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; R=$K_/rflex_${RT:?}; TAG=${TAG:?}; MD=${MSKDIAM:-220}; P=$R/post; mkdir -p $P
awk '{print NR, $3, $4, $5}' $R/leaves.txt > $P/leaves_int.txt; echo "== rflex $RT: $(grep -c . $P/leaves_int.txt) leaves"
[ -s $P/tree/merge_tree.txt ] || $PY $T/state_merge.py $P/tree --list $P/leaves_int.txt --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "states,\|^== " | head -4
$PY $T/tree_wiener4.py $P/tw --list $P/leaves_int.txt --labels $R/final_labels.txt --tree $P/tree/merge_tree.txt --gt $T/sim$TAG/gt_particles.txt --truth $T/sim$TAG/maps --mskdiam $MD 2>&1 | grep -v Warning
echo DONE
