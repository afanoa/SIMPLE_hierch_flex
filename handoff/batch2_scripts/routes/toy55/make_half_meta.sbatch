#!/bin/bash
#SBATCH -J half_meta
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan; O=$A/rnasep_klem/f692e_halfmeta; mkdir -p $O
for h in A B; do $PY $W/make_half_meta.py $A/rnasep_coupled/basis_span_all/meta_all.npz $W/simrnasep692_$h/meta_rows.txt $O/meta_half_$h.npz || exit 1; done
