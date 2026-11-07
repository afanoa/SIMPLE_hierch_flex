#!/bin/bash
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# re-render the m33 and base80 scorecard pages (scorecard_page.py only) so they open with the closed-state gate line
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for r in full692_m33_s1 full692_base80_s1; do $PY $T/scorecard_page.py $K/scorecard/$r $K/stree_port/$r/1_hetero_tree $r rnasep | head -4; done
