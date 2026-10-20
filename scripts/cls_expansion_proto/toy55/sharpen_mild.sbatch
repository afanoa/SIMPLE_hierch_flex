#!/bin/bash
#SBATCH -J rnasep_sharpen_mild
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; P=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/matched_rec
$PY $W/sharpen_matched.py $P/maps $P/for_download_mild --rmask 116 --bfac -150 | head -5
