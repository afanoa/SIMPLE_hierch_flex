#!/bin/bash
#SBATCH -J klem_p1hold2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# 1-D capture ratio and per-class gains from the p1hold scoring tables (p1_holdout2.py)
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
$PY $T/p1_holdout2.py $K_/p1hold/r22 $K_/p1hold/r11 $K_/p1hold/r12 $K_/p1hold/r21 --nb 8 --gt $T/simswing30/gt_particles.txt 2>&1 | grep -v Warning; echo DONE
