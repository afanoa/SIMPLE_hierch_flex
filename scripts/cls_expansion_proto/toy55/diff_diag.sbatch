#!/bin/bash
#SBATCH -J klem_diffdiag
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem
$PY $W/diff_diag.py $K/synth_${MODE:?}_A $K/synth_${MODE}_B ${T:-splitparent_k6} ${SA:?} ${SB:?} $W/maps160/toy_consensus.mrc $K/diffdiag_${MODE}.png 2>&1 | grep -v Warning
