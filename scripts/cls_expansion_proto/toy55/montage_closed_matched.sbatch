#!/bin/bash
#SBATCH -J rnasep_montage_closed
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; M=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/matched_rec/maps; D=$W/maps160/docked; F=$W/figs_now
$PY $W/montage_maps.py $F/closed_matched_lp55.png --lp 5.5 "$W/maps160/toy_consensus.mrc:consensus" "$M/recvol_state03.mrc:ours 18 (5.3 A)" "$(ls $D/conf18_*.mrc | head -1):deposited 18 (2.9 A)" "$(ls $D/conf27_*.mrc | head -1):deposited 27 (open)" 2>&1 | grep -v Warning
$PY $W/montage_maps.py $F/closed_matched_lp8.png --lp 8 "$W/maps160/toy_consensus.mrc:consensus" "$M/recvol_state03.mrc:ours 18" "$(ls $D/conf18_*.mrc | head -1):deposited 18" "$(ls $D/conf27_*.mrc | head -1):deposited 27 (open)" 2>&1 | grep -v Warning
