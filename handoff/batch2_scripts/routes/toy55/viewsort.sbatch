#!/bin/bash
#SBATCH -J klem_viewsort
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; for m in ${MODES:-null_200k two_50k 55_200k}; do echo "== synth_$m"; $PY $W/view_sort_check.py $A/rnasep_klem/synth_$m $A/rnasep_klem/synth_$m/meta_toy.npz 2>&1 | grep -v Warning | tail -14; done; echo "VIEWSORT DONE"
