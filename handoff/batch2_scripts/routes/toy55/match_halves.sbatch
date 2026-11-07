#!/bin/bash
#SBATCH -J match_halves
#SBATCH -p unlimited
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# env RA RB (half runs), REF (reference run, f692e), TAG (consolidated_12_<TAG>)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for t in ${TAGS:-fullw_region fullw_zoom}; do echo "== $t"; $PY $W/match_halves.py $K/${RA:?}/consolidated_12_$t $K/${RB:?}/consolidated_12_$t $K/${REF:-f692e}/consolidated_12_$t 2>&1 | grep -v Warning; done
