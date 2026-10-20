#!/bin/bash
#SBATCH -J klem_umap_real_matched
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; F=$W/figs_now
$PY $W/klem_umap.py $P/klem_fb $P/meta_fb.npz --labels $F/matched_conformer_labels.txt --other $A/rnasep_klem/expand1/fam1 --label_title 'cells matched to the deposited conformers (particle = its argmax cell)' --label_name conformer --out $F/umap_real_fam1_fb_matched.png --out2 $F/umap_real_fam1_flex_vs_klem.png 2>&1 | grep -v Warning
