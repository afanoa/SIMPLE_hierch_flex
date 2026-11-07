#!/bin/bash
#SBATCH -J klem_closed_latent2
#SBATCH -p norm
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# separability of the PUREST earlier closed set (f692r node 55 cell 25, closed-arm 1.03, 2532 particles) in r22's whole-molecule and arm-focused flex latents
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/rflex_f692rf
REF=$K_/f692r/coordtree_fam1/nodes2/node055/even_bins.txt
for lat in "whole-molecule:$(ls -td $R/node_r22/*_flex_pca | head -1)" "arm-focused:$(ls -td $K_/focused_flex_f692rf_r22/*_flex_pca | head -1)"; do echo "== ${lat%%:*} latent, reference = f692r node 55 cell 25"
  $PY $T/closed_latent.py --flex ${lat#*:}/flex_pca_coordinates.txt --oris $K_/cheapkids_f692rf_r22/ptcl3D_parent.txt --rows $R/rows_r22.txt --ref $REF --ids 25 2>&1 | grep -v Warning | grep -v "  IC"; done
echo DONE
