#!/bin/bash
#SBATCH -J rnasep_sharpen_top30k
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# B-sharpened (-150) unmasked maps of the family-1-wide docked ranking (top 30k / next 20k / bottom 30k) for download.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; FD=${FD:-$A/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal}
$PY $W/sharpen_matched.py $FD/lever/maps $FD/lever/for_download --rmask 0 --bfac -150 2>&1 | grep -v Warning | tail -4; ls -la $FD/lever/for_download | awk '{print $5, $9}'
