#!/bin/bash
#SBATCH -J klem_count_states
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Distinct reproducible states across the leaves: toy (level-3 grandchildren of the balanced realigned tree, truth
# calibrates the similarity threshold) and real family 1 (feedback board + its two level-3 grandchildren, deposited
# matches where available).
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan
echo "== TOY: level-3 leaves (4 grandchildren) vs 10 true conformers"; $PY $W/count_states.py $W/lat/klem_l3b_fam1_1 $W/lat/klem_l3b_fam1_2 $W/lat/klem_l3b_fam2_1 $W/lat/klem_l3b_fam2_2 --gt $W/sim/gt_particles.txt 2>&1 | grep -v Warning
echo "== TOY: level-2 children + level-3 leaves"; $PY $W/count_states.py $W/lat/klem_expand_realign1_bal05_fam1 $W/lat/klem_expand_realign1_bal05_fam2 $W/lat/klem_l3b_fam1_1 $W/lat/klem_l3b_fam1_2 $W/lat/klem_l3b_fam2_1 $W/lat/klem_l3b_fam2_2 --gt $W/sim/gt_particles.txt 2>&1 | grep -v Warning
echo "== REAL family 1: feedback board + level-3 grandchildren"; $PY $W/count_states.py $A/rnasep_klem/expand1/fam1_pipeline/klem_fb $A/rnasep_klem/expand1/fam1_l3/sub1 $A/rnasep_klem/expand1/fam1_l3/sub2 2>&1 | grep -v Warning
