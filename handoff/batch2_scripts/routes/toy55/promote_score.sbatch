#!/bin/bash
#SBATCH -J klem_promote
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# promoted IgM site posteriors scored like the kept sites (promote_score.py); env ENTRIES 'name:post_file:angle ..', OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
$PY $T/promote_score.py $K/factorial/armsigm1_truth/arms_truth.txt ${ENTRIES:?} 2>&1 | grep -v Warning | tee ${OUT:?}
