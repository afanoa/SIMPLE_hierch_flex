#!/bin/bash
#SBATCH -J global_embedding
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 05:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# One embedding of the whole tree for a finished run (waits for the family DONE markers). env: RTAG (f692), NSUB (80000), META (meta_all), GT (toy truth file, optional), HOLO0 (101)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692}
n=0; while [ $n -lt 240 ] && { [ ! -f $R/fam1/DONE ] || [ ! -f $R/fam2/DONE ]; }; do sleep 60; n=$((n+1)); done
$PY $W/global_embedding.py $R ${META:-$A/rnasep_coupled/basis_span_all/meta_all.npz} $R/fig_global_embedding.png --nsub ${NSUB:-80000} ${GT:+--gt $GT --holo0 ${HOLO0:-101}} 2>&1 | grep -v Warning
echo "EMBEDDING DONE"
