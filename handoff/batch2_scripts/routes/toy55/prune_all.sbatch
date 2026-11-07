#!/bin/bash
#SBATCH -J f692e_prune
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# pruning test of the f692e landscape: every board of LIST (one dir per line) with index % NCHUNK == CHUNK -> prune_one.sh
W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
i=0; while read B; do [ -n "$B" ] || continue; if [ $((i % ${NCHUNK:-1})) = ${CHUNK:-0} ]; then n=$(echo $B | sed "s#$K/f692e/##; s#/#_#g"); bash $W/prune_one.sh $B $n $K/f692e_prune/$n; fi; i=$((i + 1)); done < ${LIST:?}
echo "CHUNK ${CHUNK:-0} DONE"
