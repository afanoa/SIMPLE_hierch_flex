#!/bin/bash
#SBATCH -J klem_closed_latent
#SBATCH -p norm
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/rflex_f692rf
for set in "16,15" "16,15,9,25"; do echo "== reference f692r family-1 states $set (closed-arm 0.87 / 0.83 / 0.60 / 0.59)"
  $PY $T/closed_latent.py --flex $(ls -td $R/node_r22/*_flex_pca | head -1)/flex_pca_coordinates.txt --oris $K_/cheapkids_f692rf_r22/ptcl3D_parent.txt --rows $R/rows_r22.txt --ref $K_/f692r/coordtree_fam1/units.txt --ids $set 2>&1 | grep -v Warning; done
echo DONE
