#!/bin/bash
#SBATCH -J rnasep_montage_docked
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; D=$W/maps160/docked; RK=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank/maps
$PY $W/montage_maps.py $W/figs_now/docked_vs_undocked_lp8.png --lp 8 "$RK/recvol_state01.mrc:docked top 15k (tip 0.75)" "$RK/recvol_state03.mrc:undocked bottom 15k (tip 0.19)" "$(ls $D/conf18_*.mrc | head -1):deposited 18 (closed)" "$(ls $D/conf27_*.mrc | head -1):deposited 27 (open)" "$W/maps160/toy_consensus.mrc:consensus" 2>&1 | grep -v Warning
