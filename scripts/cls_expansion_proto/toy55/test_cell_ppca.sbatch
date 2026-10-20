#!/bin/bash
#SBATCH -J toy10_cell_ppca
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# cell-weighted PPCA vs the unweighted cell PCA vs the global flex latent, on three toy10 boards (parent, a family board,
# a level-3 feedback board), truth readouts on each board's particles
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; GT=$W/sim/gt_particles.txt; FD=$(ls -td $W/lat/flex12/*_flex_pca | head -1)
STATES=1,10,18,27,28,37,46,48,49,51; AXIS=-0.246,-0.916,-0.442,1.000,0.972,-0.685,-0.489,0.446,0.499,-0.139
for b in lat/klem lat/klem_expand_realign1_bal05_m232_fam2 lat/klem_l3b_fb_fam1_1 lat/klem_l3_fam2_2; do echo "=================== $b"; $PY $W/cell_ppca.py $W/$b $W/$b/cell_ppca --npc 4 --gt $GT --states $STATES --axis=$AXIS --flexdir $FD 2>&1 | grep -v Warning; done
echo "ALL DONE"
