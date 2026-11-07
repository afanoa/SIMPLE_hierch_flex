#!/bin/bash
#SBATCH -J rnasep_montage_unmatched
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; U=/mnt/beegfs/elmlund/afan/rnasep_klem/unmatched_states/fam2/maps; D=$W/maps160/docked
$PY $W/montage_maps.py $W/figs_now/unmatched_states_lp10.png --lp 10 "$W/maps160/toy_consensus.mrc:consensus" "$U/recvol_state01.mrc:U1 27k, 8.4 A" "$U/recvol_state02.mrc:U2 18k, 9.0 A" "$U/recvol_state03.mrc:U3 18k, 9.4 A" "$U/recvol_state04.mrc:U4 13k, 10.2 A" "$U/recvol_state05.mrc:U5 9k, 10.6 A" "$U/recvol_state06.mrc:U6 9k, 9.4 A" "$U/recvol_state07.mrc:U7 5k, 9.8 A" "$(ls $D/conf37_*.mrc | head -1):deposited 37" "$(ls $D/conf27_*.mrc | head -1):deposited 27" 2>&1 | grep -v Warning
