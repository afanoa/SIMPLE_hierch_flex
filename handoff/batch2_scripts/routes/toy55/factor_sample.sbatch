#!/bin/bash
#SBATCH -J klem_factors
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# factor sampling of run RTAG (factor_sample.py; boards GLOB, default region*), then one reconstruction + readout job per factor
# (mu_rec.sbatch + pooled_analysis.sbatch + tidy_halves.sh). env RTAG, OUT (factors), GLOB, NBIN (10), BUDGET (3000), MINEND (1000), GT, CODES
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; O=$R/${OUT:-factors}; export OMP_NUM_THREADS=8
$PY $W/factor_sample.py $R $O --glob "${GLOB:-region*}" --nbin ${NBIN:-10} --budget ${BUDGET:-3000} --minend ${MINEND:-1000} $( [ -n "${GT:-}" ] && echo "--gt $GT --codes ${CODES:-}" ) 2>&1 | grep -v Warning | tee $O.log
for d in $O/fam*_factor*; do c=${OUT:-factors}/$(basename $d)
  sbatch --parsable -J klem_factorrec -c 4 --mem=32G -t 06:00:00 -o $A/rnasep_klem/klem_factorrec.%j.out --wrap "RTAG=$RTAG CDIR=$c bash $W/mu_rec.sbatch > /dev/null 2>&1; RTAG=$RTAG CDIR=$c bash $W/pooled_analysis.sbatch 2>&1 | grep -A14 'arm readout'; bash $W/tidy_halves.sh $R/$c"
done
